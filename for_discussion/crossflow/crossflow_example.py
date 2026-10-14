#!/usr/bin/env python3
"""Synthetic example: how the fracture connection factors give crossflow in Flow.

Generates the figures and the numbers (numbers.tex) used by crossflow_example.tex.
All potentials are in bar, rates in reservoir m3/day, conductances in m3/day/bar.
"""
import os
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figs")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3})

# ---------------------------------------------------------------- parameters
g = 9.80665
bar = 1e5
day = 86400.0
CONV = bar * day                 # m3/(s Pa) -> m3/(d bar)
mD = 9.869233e-16

rho_r = 1020.0                   # reservoir (formation) water density, fracture fluid in the model
mu = 0.5e-3                      # water viscosity [Pa s]
lam = 1.0 / mu                   # water mobility [1/(Pa s)]
perm = 100 * mD                  # permeability
h = 10.0                         # layer thickness
dx = 20.0                        # reservoir cell size (x = y)
d_leak = 0.25 * dx               # leak-off distance, legacy model: 0.25 * cell size
sides = 2                        # legacy leak-off model: both fracture faces
rw = 0.1                         # wellbore radius
R = 12.0                         # fracture radius
zo = 2000.0                      # seed depth = centre of middle completion
zc = np.array([1990.0, 2000.0, 2010.0])   # completion / cell centre depths
names = ["top", "middle", "bottom"]
z_ref = 1985.0                   # well reference (BHP) depth
p_o = 200.0                      # reservoir pressure at z_o [bar]
E = 10e9                         # Young's modulus
nu = 0.25
sig_grad = 0.17                  # minimum horizontal stress gradient [bar/m]
p_net0 = 2.0                     # net pressure at the seed [bar]
Q_base = 8.0                    # "fairly low" injection rate [m3/d]
gamma_base = 0.004               # depth-linear potential mismatch [bar/m]
W_over_L = 100.0                 # well -> fracture conductance relative to total leak-off


def band_area(R, a, b):
    """Area of the disk r < R between z = a and z = b (relative to the centre)."""
    a, b = max(a, -R), min(b, R)
    if b <= a:
        return 0.0
    F = lambda z: z * np.sqrt(R * R - z * z) + R * R * np.arcsin(z / R)
    return F(b) - F(a)


A = np.array([band_area(R, z - h / 2 - zo, z + h / 2 - zo) for z in zc])
Lk = sides * lam * perm * A / d_leak * CONV                  # leak-off conductance per layer
r_o = 0.198 * dx                                             # Peaceman radius, square cell
Mk = np.full(3, 2 * np.pi * perm * h * lam / np.log(r_o / rw) * CONV)  # matrix CTF*lambda
W = W_over_L * Lk.sum()


# ------------------------------------------------- 2D fracture model solution
def fracture_2d(phi_layers, phi_w, n=96):
    """Steady fracture flow on a Cartesian grid in the fracture plane (x, z).

    Unknown: potential Phi = p - rho_r g z (bar), the quantity the fracture code uses
    (fracture_dgh_ = g rho z).  Cubic-law conductivity, leak-off to the layer holding
    the cell centre, ring feed from the well with total conductance W.
    """
    dl = 2 * R / n
    xs = -R + dl * (np.arange(n) + 0.5)
    X, Zr = np.meshgrid(xs, xs)                 # Zr relative to z_o, positive down
    r = np.hypot(X, Zr)
    inside = r < R
    idx = -np.ones(X.shape, int)
    idx[inside] = np.arange(inside.sum())
    nc = inside.sum()
    wmax = 8 * (1 - nu**2) / (np.pi * E) * p_net0 * bar * R
    width = np.where(inside, 8 * (1 - nu**2) / (np.pi * E) * p_net0 * bar
                     * np.sqrt(np.clip(R**2 - r**2, 0, None)), 0.0)
    width = np.maximum(width, 1e-4)
    rows, cols, vals = [], [], []
    diag = np.zeros(nc)
    rhs = np.zeros(nc)
    # faces (cubic law, harmonic mean of half transmissibilities)
    for (di, dj) in ((0, 1), (1, 0)):
        a = inside[:n - di, :n - dj] & inside[di:, dj:]
        i1 = idx[:n - di, :n - dj][a]
        i2 = idx[di:, dj:][a]
        w1 = width[:n - di, :n - dj][a]
        w2 = width[di:, dj:][a]
        t = lam * CONV / (12 / w1**3 + 12 / w2**3) * 2  # dl/dl * face, both halves
        rows += [i1, i2]; cols += [i2, i1]; vals += [-t, -t]
        np.add.at(diag, i1, t); np.add.at(diag, i2, t)
    # leak-off
    zabs = zo + Zr[inside]
    layer = np.clip(np.floor((zabs - (zc[0] - h / 2)) / h).astype(int), 0, 2)
    lcell = sides * lam * perm * dl * dl / d_leak * CONV
    diag += lcell
    rhs += lcell * phi_layers[layer]
    # ring feed
    ring = r[inside] < max(0.6, 1.5 * dl)
    wcell = W / ring.sum()
    diag[ring] += wcell
    rhs[ring] += wcell * phi_w
    Amat = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                         shape=(nc, nc)) + sp.diags(diag)
    phi = spla.spsolve(Amat.tocsc(), rhs)
    leak = lcell * (phi - phi_layers[layer])               # m3/d per cell
    q = np.array([leak[layer == k].sum() for k in range(3)])
    phiF = np.array([np.average(phi[layer == k]) for k in range(3)])
    return dict(X=X, Z=Zr + zo, inside=inside, idx=idx, phi=phi, width=width,
                leak=leak / (dl * dl), q=q, phiF=phiF, dl=dl, wmax=wmax,
                feed=(wcell * (phi_w - phi[ring])).sum())


def field(res, values):
    out = np.full(res["X"].shape, np.nan)
    out[res["inside"]] = values
    return out



# ------------------------------------------------- linear response of the fracture model
def _response():
    """q_k and leak-weighted fracture potential per layer are linear in (Phi_w, Phi_layers)."""
    base = fracture_2d(np.zeros(3), 1.0)
    G, a = base["q"].copy(), base["phiF"].copy()
    H = np.zeros((3, 3)); B = np.zeros((3, 3))
    for j in range(3):
        e = np.zeros(3); e[j] = 1.0
        r = fracture_2d(e, 0.0)
        H[:, j], B[:, j] = r["q"], r["phiF"]
    return G, H, a, B


G_resp, H_resp, a_resp, B_resp = _response()


def fracture_model(phi_layers, phi_w):
    """Fracture solve (perf_pressure control at the seed): leak-off per layer and the
    leak-off weighted fracture potential per layer."""
    q = G_resp * phi_w + H_resp @ phi_layers
    phiF = a_resp * phi_w + B_resp @ phi_layers
    return phiF, q


def legacy_ctf(q, phi_w, phi_layers):
    """Fracture::wellIndices_ (wi_upscaling = legacy): WI = q / dp, negative -> 0."""
    dp = phi_w - phi_layers
    safe = np.where(np.abs(dp) > 1e-14, dp, 1.0)
    T = np.where(np.abs(dp) > 1e-14, q / safe, 0.0)
    return np.where(T > 0, T, 0.0), dp


P_FLOOR = 1e4 / bar     # solver.wi_pressure_floor (default 1e4 Pa) in bar
ALPHA_MAX = 2.0          # solver.wi_normalization_max


def conductivity_ctf(q, phi_w, phi_layers, norm=True, gate=False, floor=P_FLOOR):
    """Fracture::wellIndices_ (wi_upscaling = conductivity).

    ctf = sum leakof/mobility (T_k = L_k), optional sign gate, and one global factor
    alpha = min(sum q / sum T max(dp, floor), alpha_max) when wi_flux_normalization is on.
    """
    dp = phi_w - phi_layers
    T = Lk.copy()
    not_fed = (q < 0.0) | ((dp <= 0.0) & (q > 0.0))
    if gate:
        T = np.where(not_fed, 0.0, T)
    sum_q = q.sum()
    sum_ctf_dp = (T * np.maximum(dp, floor)).sum()
    alpha = 1.0
    if norm and sum_q > 0.0 and sum_ctf_dp > 0.0:
        alpha = min(sum_q / sum_ctf_dp, ALPHA_MAX)
    return alpha * T, alpha


def flow_rate_control(T, Q, phi_layers, s, off=None):
    """Flow's standard-well model, rate control, crossflow allowed.

    q_k = T_k (Phi_w + o_k - Phi_k) + M_k (Phi_w + s_k - Phi_k), with o_k = s_k for the
    code's CTFs and o_k = the fracture offset for option 3.
    """
    o = s if off is None else off
    phi_w = (Q - (T * (o - phi_layers)).sum() - (Mk * (s - phi_layers)).sum()) / (T + Mk).sum()
    return phi_w, T * (phi_w + o - phi_layers), Mk * (phi_w + s - phi_layers)


def true_state(Q, phi_layers):
    """Reference: fracture model + matrix completions at a common well potential, s = 0."""
    phi_w = (Q - (H_resp @ phi_layers).sum() + (Mk * phi_layers).sum()) / (G_resp.sum() + Mk.sum())
    phiF, qf = fracture_model(phi_layers, phi_w)
    return phi_w, phiF, qf, Mk * (phi_w - phi_layers)


def cond_mode_ctf(mode, q, phi_w, phi_layers, s):
    """cond: defaults (normalisation, floor); cond_nonorm: wi_flux_normalization=false;
    cond_nofloor: floor -> 0; cond_gate_col: sign gate + wi_flow_column."""
    if mode == "cond":
        return conductivity_ctf(q, phi_w, phi_layers)
    if mode == "cond_nonorm":
        return conductivity_ctf(q, phi_w, phi_layers, norm=False)
    if mode == "cond_nofloor":
        return conductivity_ctf(q, phi_w, phi_layers, floor=1e-12)
    if mode == "cond_gate_col":
        return conductivity_ctf(q, phi_w + s, phi_layers, gate=True)
    raise ValueError(mode)


def coupled_code(Q, phi_layers, s, mode="legacy", iters=200):
    """Fixed point of the sequential coupling (fracture solve <-> Flow well solve).

    The fracture is solved at the potential of Flow's seed connection; the CTFs are formed
    from it and handed to Flow; Flow's rate control gives a new seed potential.
    mode: legacy | flow_column (option 2) | fracture_pressure (option 3).
    """
    phi_w = Q / (Lk.sum() + Mk.sum())
    T = Lk.copy()
    for _ in range(iters):
        phiF, q = fracture_model(phi_layers, phi_w)
        if mode == "legacy":
            T, _ = legacy_ctf(q, phi_w, phi_layers)
            phi_w, qf, qm = flow_rate_control(T, Q, phi_layers, s)
        elif mode == "flow_column":
            T, _ = legacy_ctf(q, phi_w + s, phi_layers)
            phi_w, qf, qm = flow_rate_control(T, Q, phi_layers, s)
        elif mode.startswith("cond"):
            T, _ = cond_mode_ctf(mode, q, phi_w, phi_layers, s)
            phi_w, qf, qm = flow_rate_control(T, Q, phi_layers, s)
        else:
            T = Lk
            phi_w, qf, qm = flow_rate_control(T, Q, phi_layers, s, off=phiF - phi_w)
    return phi_w, qf, qm, T


# =============================================================== computations
num = {}
phi0 = np.zeros(3)                     # hydrostatic reservoir: equal layer potentials
s_base = -gamma_base * (zc - zo) + 0.0  # mismatch: +top, -bottom

# reference state at Q_base
phiw_ref, phiF_ref_k, qf_ref, qm_ref = true_state(Q_base, phi0)
phiF_ref = float(np.average(phiF_ref_k, weights=Lk))
res2d = fracture_2d(phi0, phiw_ref)
T_leg, dp_leg = legacy_ctf(res2d["q"], phiw_ref, phi0)

# Flow with these CTFs, no mismatch / with mismatch (same CTFs)
phiw_A, qfA, qmA = flow_rate_control(T_leg, Q_base, phi0, np.zeros(3))
phiw_C, qfC, qmC = flow_rate_control(T_leg, Q_base, phi0, s_base)
Q_crit = ((T_leg + Mk) * (s_base - s_base.min())).sum()

# remedies at Q_base (fixed points of the sequential coupling)
rem = {m: coupled_code(Q_base, phi0, s_base, m) for m in ("legacy", "flow_column", "fracture_pressure")}

# connection at seed depth (before the depth fix): s_k = -rho g (z_k - z_o)
s_seed = -rho_r * g * (zc - zo) / bar
Q_crit_seed = ((T_leg + Mk) * (s_seed - s_seed.min())).sum()
_, qf_seed, qm_seed = flow_rate_control(T_leg, Q_base, phi0, s_seed)

# Q sweep
Qs = np.linspace(0.5, 40, 160)
true_q = np.array([sum(true_state(Q, phi0)[2:]) for Q in Qs])
code_q = np.array([sum(coupled_code(Q, phi0, s_base)[1:3]) for Q in Qs])
code_Qcrit = Qs[np.argmax(code_q[:, 2] >= 0)]

# gamma sweep for Q_crit
gammas = np.linspace(0, 0.02, 81)
Qc_g = [((T_leg + Mk) * (-(gm) * (zc - zo) + gm * (zc[-1] - zo))).sum() for gm in gammas]

# physical non-hydrostatic layers (scenario B): legacy CTF singularity
delta_b = 0.03 * np.array([-1.0, 0.0, 1.0])     # bottom layer at higher potential
phiws = np.linspace(-0.02, 0.12, 2801)
Tb = []
for pw in phiws:
    _, qq = fracture_model(delta_b, pw)
    dp = pw - delta_b
    Tb.append(qq / np.where(np.abs(dp) > 1e-9, dp, np.nan))
Tb = np.array(Tb)
# q_b = 0 at G_b pw + H_b . delta = 0
pw_qb0 = -(H_resp[2] @ delta_b) / G_resp[2]

# dynamic sequential coupling with reservoir storage (tank per layer)
C = np.full(3, 8.0)            # m3/bar, layer storage seen by the connection
Aq = np.full(3, 1000.0)         # m3/d/bar, pressure support from the rest of the layer
dt = 1.0
nstep = 40


def dynamic(mode, s):
    phi = phi0.copy()
    phi_w = Q_base / (Lk.sum() + Mk.sum())
    T_old = None
    hist = []
    for n in range(nstep):
        if mode == "true":
            nn = 4
            Am = np.zeros((nn, nn)); b = np.zeros(nn)
            for i in range(3):
                # C(phi-phi0)/dt + Aq phi = q_i + M_i(phi_w - phi_i), q = G phi_w + H phi
                Am[i, :3] = -H_resp[i]
                Am[i, i] += C[i] / dt + Aq[i] + Mk[i]
                Am[i, 3] = -(G_resp[i] + Mk[i])
                b[i] = C[i] / dt * phi[i]
            Am[3, :3] = H_resp.sum(axis=0) - Mk
            Am[3, 3] = G_resp.sum() + Mk.sum()
            b[3] = Q_base
            x = np.linalg.solve(Am, b)
            phi, phi_w = x[:3], x[3]
            hist.append(G_resp * phi_w + H_resp @ phi + Mk * (phi_w - phi))
            continue
        phiF, q = fracture_model(phi, phi_w)
        if mode == "fracture_pressure":
            Tf = Lk
            of = phiF - phi_w
        elif mode.startswith("cond"):
            Tn, _ = cond_mode_ctf(mode, q, phi_w, phi, s)
            Tf = Tn if T_old is None else (Tn + 2.0 * T_old) / 3.0   # wellIndicesAvrg
            T_old = Tf
            of = s
        else:
            dps = (phi_w + s) if mode == "flow_column" else phi_w
            Tn, _ = legacy_ctf(q, dps, phi)
            Tf = Tn if T_old is None else (Tn + 2.0 * T_old) / 3.0   # wellIndicesAvrg
            T_old = Tf
            of = s
        nn = 4
        Am = np.zeros((nn, nn)); b = np.zeros(nn)
        for i in range(3):
            Am[i, i] = C[i] / dt + Aq[i] + Tf[i] + Mk[i]
            Am[i, 3] = -(Tf[i] + Mk[i])
            b[i] = C[i] / dt * phi[i] + Tf[i] * of[i] + Mk[i] * s[i]
        Am[3, 3] = (Tf + Mk).sum(); Am[3, :3] = -(Tf + Mk)
        b[3] = Q_base - (Tf * of).sum() - (Mk * s).sum()
        x = np.linalg.solve(Am, b)
        phi, phi_w = x[:3], x[3]
        hist.append(Tf * (phi_w + of - phi) + Mk * (phi_w + s - phi))
    return np.array(hist)


dyn = {m: dynamic(m, s_base) for m in ("true", "legacy", "flow_column", "fracture_pressure")}

# ---------------------------------------------------- conductivity mode
cond_modes = ("cond_nonorm", "cond", "cond_nofloor", "cond_gate_col")
cond_fp = {m: coupled_code(Q_base, phi0, s_base, m) for m in cond_modes}
cond_alpha = {}
cond_fracshare = {}
for m in cond_modes:
    pw, qf_, qm_, T_ = cond_fp[m]
    _, qmod = fracture_model(phi0, pw)
    cond_alpha[m] = (T_ / Lk).max()
    cond_fracshare[m] = qf_.sum() / qmod.sum()     # Flow's fracture injection / model leak-off
cond_q = {m: np.array([sum(coupled_code(Q, phi0, s_base, m)[1:3]) for Q in Qs]) for m in cond_modes}
cond_Qcrit = {m: Qs[np.argmax(cond_q[m][:, 2] >= 0)] for m in cond_modes}
# analytic: no normalisation -> T = L
Q_crit_nonorm = ((Lk + Mk) * (s_base - s_base.min())).sum()
# fixed point with floor: alpha = g_s phi_w / (sum L floor), phi_w (alpha sum L + sum M) = Q
g_s = G_resp.sum()
# crossflow onset at phi_w = -s_min
pw_on = -s_base.min()
Q_crit_floor = g_s * pw_on**2 / P_FLOOR + Mk.sum() * pw_on
# alpha vs phi_w (hydrostatic layers, no mismatch)
alpha_curve = np.array([conductivity_ctf(fracture_model(phi0, pw)[1], pw, phi0)[1] for pw in phiws[phiws > 0]])
alpha_nofloor = np.array([conductivity_ctf(fracture_model(phi0, pw)[1], pw, phi0, floor=1e-12)[1]
                          for pw in phiws[phiws > 0]])
share_Q = np.array([coupled_code(Q, phi0, np.zeros(3), "cond")[1].sum()
                    / fracture_model(phi0, coupled_code(Q, phi0, np.zeros(3), "cond")[0])[1].sum()
                    for Q in Qs])
# scenario B with conductivity: alpha and gate
condB = []
for pw in phiws:
    _, qq = fracture_model(delta_b, pw)
    Tc, al = conductivity_ctf(qq, pw, delta_b)
    Tg, _ = conductivity_ctf(qq, pw, delta_b, gate=True)
    condB.append(np.concatenate([Tc / Lk, Tg / Lk]))
condB = np.array(condB)
dyn_cond = {m: dynamic(m, s_base) for m in ("cond", "cond_nonorm", "cond_gate_col")}

# ================================================================== figures
cols = ["tab:blue", "tab:green", "tab:red"]

# Fig 1: geometry
fig, ax = plt.subplots(figsize=(6.2, 4.2))
for k in range(3):
    ax.add_patch(Rectangle((-20, zc[k] - h / 2), 40, h, fc=["#f3e9d2", "#efe0bd", "#e8d5a6"][k],
                           ec="k", lw=0.5))
    ax.text(-19.3, zc[k], f"layer {k+1} ({names[k]})\n$z_{k+1}$ = {zc[k]:.0f} m",
            va="center", fontsize=8)
ax.add_patch(Circle((0, zo), R, fc="tab:orange", alpha=0.35, ec="tab:orange", lw=1.5))
ax.plot([0, 0], [1978, 2016], color="k", lw=3)
for k in range(3):
    ax.plot([-0.9, 0.9], [zc[k]] * 2, color="tab:red", lw=6, solid_capstyle="butt")
    ax.annotate(f"completion {k+1}\n$A_{k+1}$ = {A[k]:.0f} m$^2$", (0.9, zc[k]),
                (13.5, zc[k]), fontsize=8, va="center",
                arrowprops=dict(arrowstyle="-", lw=0.5))
ax.plot(0, zo, "k*", ms=10)
ax.annotate("seed $z_o$", (0, zo), (4, zo - 3.2), fontsize=8)
ax.annotate("", (R, zo + 0.0), (0, zo), arrowprops=dict(arrowstyle="<->", lw=0.8))
ax.text(R / 2, zo + 1.3, f"R = {R:.0f} m", ha="center", fontsize=8)
ax.text(0.6, 1979.5, "vertical injector", fontsize=8)
ax.set_xlim(-20, 25); ax.set_ylim(2017, 1978)
ax.set_xlabel("horizontal distance in the fracture plane [m]"); ax.set_ylabel("depth z [m]")
ax.set_aspect("equal"); ax.grid(False)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "geometry.pdf")); plt.close(fig)

# Fig 2: fracture model fields
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.7))
Z = res2d["Z"]; X = res2d["X"]
p_field = field(res2d, res2d["phi"]) + p_o + rho_r * g * (Z - zo) / bar
items = [(p_field, "fracture pressure $p_f$ [bar]", "viridis"),
         (field(res2d, res2d["width"][res2d["inside"]]) * 1e3, "width $w$ [mm]", "magma"),
         (field(res2d, res2d["leak"]) * 1e3, "leak-off flux [L/(d m$^2$)]", "cividis")]
for ax, (v, title, cm) in zip(axs, items):
    pc = ax.pcolormesh(X, Z, v, shading="auto", cmap=cm, edgecolors="none", linewidth=0, rasterized=True)
    for zb in (zc[0] + h / 2, zc[1] + h / 2):
        ax.axhline(zb, color="w", lw=0.6, ls="--")
    ax.set_ylim(zo + R, zo - R); ax.set_aspect("equal"); ax.set_title(title, fontsize=8)
    ax.set_xlabel("x [m]"); ax.grid(False)
    fig.colorbar(pc, ax=ax, shrink=0.8)
axs[0].set_ylabel("depth z [m]")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fracture_fields.pdf")); plt.close(fig)

# Fig 3: potentials vs depth
zz = np.linspace(zo - R, zo + R, 200)
fig, ax = plt.subplots(figsize=(5.6, 3.4))
for k in range(3):
    ax.axhspan(zc[k] - h / 2, zc[k] + h / 2, color=["#f3e9d2", "#efe0bd", "#e8d5a6"][k], zorder=0)
ax.axvline(0, color="k", lw=1.0, label=r"reservoir layers $\Phi_k$ (hydrostatic)")
ax.plot(np.full_like(zz, phiw_ref), zz, color="tab:gray", ls="--",
        label=r"well at the seed $\Phi_w$ (code: $p_{inj}-\rho_r g z_o$)")
ax.plot(field(res2d, res2d["phi"])[:, res2d["X"].shape[1] // 2], res2d["Z"][:, 0],
        color="tab:orange", lw=1.5, label=r"fracture $\Phi_f(z)$ along the centre line")
ax.plot(phiw_C + s_base[1] - gamma_base * (zz - zo), zz, color="tab:red",
        label=r"Flow's connection potential $\Phi_w + s(z)$")
for k in range(3):
    ax.plot(phiw_C + s_base[k], zc[k], "o", color="tab:red")
    ax.annotate(f"$s_{k+1}$ = {s_base[k]:+.3f}", (phiw_C + s_base[k], zc[k]),
                (phiw_C + s_base[k] + 0.006, zc[k] + 1.6), fontsize=7, color="tab:red")
ax.fill_betweenx(zz, 0, np.minimum(phiw_C - gamma_base * (zz - zo), 0), color="tab:red", alpha=0.25)
ax.set_ylim(zo + R, zo - R)
ax.set_xlabel(r"potential $\Phi = p - p_o - \rho_r g (z - z_o)$ [bar]")
ax.set_ylabel("depth z [m]")
ax.legend(fontsize=7, loc="lower right")
ax.set_title(r"$Q$ = %g m$^3$/d, $\gamma$ = %g bar/m" % (Q_base, gamma_base), fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "profiles.pdf")); plt.close(fig)

# Fig 4: rates vs Q
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True)
for ax, data, title in ((axs[0], true_q, "fracture model (reference)"),
                        (axs[1], code_q, "Flow with code CTFs, mismatch $s_k$")):
    for k in range(3):
        ax.plot(Qs, data[:, k], color=cols[k], ls=["-", "-", "--"][k],
                label=f"completion {k+1} ({names[k]})")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel(r"well rate $Q$ [m$^3$/d]"); ax.set_title(title, fontsize=8)
axs[1].axvline(Q_crit, color="k", ls=":")
axs[1].text(Q_crit + 0.5, axs[1].get_ylim()[1] * 0.85, r"$Q_{crit}$", fontsize=8)
axs[1].axvspan(0, Q_crit, color="tab:red", alpha=0.08)
axs[0].set_ylabel(r"connection rate [m$^3$/d]  (+ injection)"); axs[0].legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "rates_vs_Q.pdf")); plt.close(fig)

# Fig 5: Q_crit vs gamma
fig, ax = plt.subplots(figsize=(4.2, 2.8))
ax.plot(gammas * 100, Qc_g, "k")
ax.axvline(gamma_base * 100, color="tab:red", ls=":")
ax.set_xlabel(r"mismatch gradient $\gamma$ [bar/100 m]"); ax.set_ylabel(r"$Q_{crit}$ [m$^3$/d]")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "qcrit.pdf")); plt.close(fig)

# Fig 6: legacy CTF singularity (scenario B)
fig, ax = plt.subplots(figsize=(5.2, 3.0))
ax.plot(phiws, Tb[:, 2] / Lk[2], color="tab:red", label="bottom: raw $q_3/\\Delta p_3$")
ax.plot(phiws, np.where(Tb[:, 2] > 0, Tb[:, 2], 0) / Lk[2], color="k", lw=1.2, ls="--",
        label="bottom: legacy (negative $\\to$ 0)")
ax.plot(phiws, Tb[:, 0] / Lk[0], color="tab:blue", label="top")
ax.axvline(delta_b[2], color="gray", ls=":")
ax.axvline(pw_qb0, color="gray", ls="--")
ax.axvspan(delta_b[2], pw_qb0, color="tab:red", alpha=0.08)
ax.set_ylim(-6, 8)
ax.set_xlabel(r"well potential $\Phi_w$ [bar] (layer potentials $-0.03, 0, +0.03$)")
ax.set_ylabel(r"$T_k / L_k$")
ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "legacy_singularity.pdf")); plt.close(fig)

# Fig 7: dynamics
fig, axs = plt.subplots(1, 4, figsize=(7.4, 2.6), sharey=True)
titles = {"true": "reference", "legacy": "legacy CTF", "flow_column": "option 2",
          "fracture_pressure": "option 3"}
for ax, m in zip(axs, ("true", "legacy", "flow_column", "fracture_pressure")):
    for k in range(3):
        ax.plot(np.arange(1, nstep + 1), dyn[m][:, k], color=cols[k], label=names[k])
    ax.axhline(0, color="k", lw=0.6); ax.set_title(titles[m], fontsize=8)
    ax.set_xlabel("time [d]")
axs[0].set_ylabel(r"connection rate [m$^3$/d]"); axs[0].legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "dynamics.pdf")); plt.close(fig)

# Fig 8: remedies at Q_base (bar chart)
fig, ax = plt.subplots(figsize=(5.6, 2.8))
labels = ["fracture model", "legacy", "option 2\n(wi_flow_column)", "option 3\n(wi_fracture_pressure)"]
data = [qf_ref + qm_ref] + [rem[m][1] + rem[m][2] for m in ("legacy", "flow_column", "fracture_pressure")]
xb = np.arange(len(labels))
for k in range(3):
    ax.bar(xb + (k - 1) * 0.25, [dd[k] for dd in data], 0.25, color=cols[k], label=names[k])
ax.axhline(0, color="k", lw=0.6)
ax.set_xticks(xb); ax.set_xticklabels(labels, fontsize=7)
ax.set_ylabel(r"connection rate [m$^3$/d]"); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "remedies.pdf")); plt.close(fig)

# Fig 9: conductivity rates vs Q
fig, axs = plt.subplots(1, 4, figsize=(7.6, 2.8), sharey=True)
panels = (("legacy", code_q, Q_crit), ("conductivity,\nno normalisation", cond_q["cond_nonorm"], cond_Qcrit["cond_nonorm"]),
          ("conductivity\n(defaults)", cond_q["cond"], cond_Qcrit["cond"]),
          ("conductivity + sign gate\n+ wi_flow_column", cond_q["cond_gate_col"], cond_Qcrit["cond_gate_col"]))
for ax, (title, data, qc) in zip(axs, panels):
    for k in range(3):
        ax.plot(Qs, data[:, k], color=cols[k], ls=["-", "-", "--"][k], label=names[k])
    ax.axhline(0, color="k", lw=0.6)
    if data[0, 2] < 0:
        ax.axvspan(0, qc, color="tab:red", alpha=0.08)
    ax.set_title(title, fontsize=7.5); ax.set_xlabel(r"$Q$ [m$^3$/d]")
axs[0].set_ylabel(r"connection rate [m$^3$/d]"); axs[0].legend(fontsize=7)
axs[0].set_ylim(-6, 14)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "cond_rates_vs_Q.pdf")); plt.close(fig)

# Fig 10: alpha and fracture share
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.8))
ax = axs[0]
ax.plot(phiws[phiws > 0], alpha_curve, "k", label="floor 0.1 bar (default)")
ax.plot(phiws[phiws > 0], alpha_nofloor, "k--", label="no floor")
ax.axvline(P_FLOOR, color="gray", ls=":")
ax.set_xlabel(r"$\Phi_w$ [bar]"); ax.set_ylabel(r"$\alpha$"); ax.legend(fontsize=7)
ax.set_title(r"normalisation factor $\alpha$ (hydrostatic layers)", fontsize=8)
ax = axs[1]
ax.plot(Qs, share_Q, "k")
ax.set_xlabel(r"$Q$ [m$^3$/d]"); ax.set_ylabel("Flow fracture injection / model leak-off")
ax.set_title("consistency of the fracture flux (no mismatch)", fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "cond_alpha.pdf")); plt.close(fig)

# Fig 11: scenario B with conductivity
fig, ax = plt.subplots(figsize=(5.2, 3.0))
ax.plot(phiws, np.where(Tb[:, 2] > 0, Tb[:, 2], 0) / Lk[2], color="k", lw=1.0, ls="--",
        label="legacy, bottom")
ax.plot(phiws, condB[:, 2], color="tab:red", label="conductivity, bottom")
ax.plot(phiws, condB[:, 5], color="tab:red", ls=":", lw=1.6, label="conductivity + gate, bottom")
ax.plot(phiws, condB[:, 0], color="tab:blue", label="conductivity, top")
ax.axvline(delta_b[2], color="gray", ls=":"); ax.axvline(pw_qb0, color="gray", ls="--")
ax.set_ylim(-0.2, 3); ax.set_xlabel(r"well potential $\Phi_w$ [bar]"); ax.set_ylabel(r"$T_k/L_k$")
ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "cond_singularity.pdf")); plt.close(fig)

# Fig 12: conductivity dynamics
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.6), sharey=True)
for ax, m, title in zip(axs, ("cond_nonorm", "cond", "cond_gate_col"),
                        ("no normalisation", "defaults", "gate + wi_flow_column")):
    for k in range(3):
        ax.plot(np.arange(1, nstep + 1), dyn_cond[m][:, k] if m in dyn_cond else dynamic(m, s_base)[:, k],
                color=cols[k], label=names[k])
    ax.axhline(0, color="k", lw=0.6); ax.set_title("conductivity: " + title, fontsize=8)
    ax.set_xlabel("time [d]")
axs[0].set_ylabel(r"connection rate [m$^3$/d]"); axs[0].legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "cond_dynamics.pdf")); plt.close(fig)

# ================================================================== numbers
def f(x, nd=2):
    return f"{x:.{nd}f}"


macros = {
    "Rfrac": f(R, 0), "hlay": f(h, 0), "dxcell": f(dx, 0), "dleak": f(d_leak, 0),
    "permmD": "100", "muCp": "0.5", "rhor": f(rho_r, 0), "Qbase": f(Q_base, 0),
    "gammabase": f(gamma_base, 3), "gammabaseh": f(gamma_base * 100, 1),
    "WoverL": f(W_over_L, 0), "Emod": "10", "nupois": f(nu, 2), "pnet": f(p_net0, 1),
    "wmax": f(res2d["wmax"] * 1e3, 2), "sigmagrad": f(sig_grad, 2),
    "Atop": f(A[0], 1), "Amid": f(A[1], 1), "Abot": f(A[2], 1), "Atot": f(A.sum(), 1),
    "Ltop": f(Lk[0], 1), "Lmid": f(Lk[1], 1), "Lbot": f(Lk[2], 1), "Ltot": f(Lk.sum(), 1),
    "Mk": f(Mk[0], 1), "ro": f(r_o, 2), "Wcond": f(W, 0),
    "phiwref": f(phiw_ref, 4), "phiFref": f(phiF_ref, 4),
    "qftop": f(res2d["q"][0]), "qfmid": f(res2d["q"][1]), "qfbot": f(res2d["q"][2]),
    "qfnodetop": f(qf_ref[0]), "qfnodemid": f(qf_ref[1]), "qfnodebot": f(qf_ref[2]),
    "qmtop": f(qm_ref[0]), "qmmid": f(qm_ref[1]), "qmbot": f(qm_ref[2]),
    "Ttop": f(T_leg[0], 1), "Tmid": f(T_leg[1], 1), "Tbot": f(T_leg[2], 1),
    "TLtop": f(T_leg[0] / Lk[0], 3), "TLmid": f(T_leg[1] / Lk[1], 3), "TLbot": f(T_leg[2] / Lk[2], 3),
    "dptop": f(dp_leg[0], 4), "dpmid": f(dp_leg[1], 4), "dpbot": f(dp_leg[2], 4),
    "phiwA": f(phiw_A, 4), "qAtop": f(qfA[0] + qmA[0]), "qAmid": f(qfA[1] + qmA[1]),
    "qAbot": f(qfA[2] + qmA[2]),
    "sTopv": f(s_base[0], 3), "sBotv": f(s_base[2], 3),
    "phiwC": f(phiw_C, 4),
    "qCftop": f(qfC[0]), "qCfmid": f(qfC[1]), "qCfbot": f(qfC[2]),
    "qCmtop": f(qmC[0]), "qCmmid": f(qmC[1]), "qCmbot": f(qmC[2]),
    "qCtop": f(qfC[0] + qmC[0]), "qCmid": f(qfC[1] + qmC[1]), "qCbot": f(qfC[2] + qmC[2]),
    "Qcrit": f(Q_crit, 1), "Ttot": f((T_leg + Mk).sum(), 1),
    "sseedtop": f(s_seed[0], 3), "Qcritseed": f(Q_crit_seed, 0),
    "qseedtop": f(qf_seed[0] + qm_seed[0], 1), "qseedmid": f(qf_seed[1] + qm_seed[1], 1),
    "qseedbot": f(qf_seed[2] + qm_seed[2], 1),
    "pwqbzero": f(pw_qb0, 4), "deltab": f(delta_b[2], 2),
    "remLtop": f(sum(rem["legacy"][1:3])[0]), "remLbot": f(sum(rem["legacy"][1:3])[2]),
    "remTwotop": f(sum(rem["flow_column"][1:3])[0]), "remTwomid": f(sum(rem["flow_column"][1:3])[1]),
    "remTwobot": f(sum(rem["flow_column"][1:3])[2]),
    "remTwoTbot": f(rem["flow_column"][3][2], 1),
    "remThrtop": f(sum(rem["fracture_pressure"][1:3])[0]), "remThrmid": f(sum(rem["fracture_pressure"][1:3])[1]),
    "remThrbot": f(sum(rem["fracture_pressure"][1:3])[2]),
    "remThrfbot": f(rem["fracture_pressure"][1][2]), "remThrmbot": f(rem["fracture_pressure"][2][2]),
    "dynLbot": f(dyn["legacy"][-1, 2]), "dynLtop": f(dyn["legacy"][-1, 0]),
    "dynTbot": f(dyn["true"][-1, 2]),
    "Cstore": f(C[0], 0), "Aqsup": f(Aq[0], 0), "codeQcrit": f(code_Qcrit, 1),
    "dpdrop": f(phiw_ref - phiF_ref, 4), "dynFbot": f(dyn["fracture_pressure"][-1, 2]),
    "dynTwobot": f(dyn["flow_column"][-1, 2]),
}
def rates(m):
    v = cond_fp[m]
    return v[1] + v[2]


for m, tag in (("cond_nonorm", "CN"), ("cond", "CD"), ("cond_nofloor", "CF"), ("cond_gate_col", "CG")):
    r = rates(m)
    macros[f"q{tag}top"] = f(r[0]); macros[f"q{tag}mid"] = f(r[1]); macros[f"q{tag}bot"] = f(r[2])
    macros[f"alpha{tag}"] = f(cond_alpha[m], 3)
    macros[f"phiw{tag}"] = f(cond_fp[m][0], 4)
    macros[f"share{tag}"] = f(cond_fracshare[m], 2)
    macros[f"Qcrit{tag}"] = f(cond_Qcrit[m], 1)
    macros[f"qf{tag}bot"] = f(cond_fp[m][1][2]); macros[f"qm{tag}bot"] = f(cond_fp[m][2][2])
macros["QcritNonormAn"] = f(Q_crit_nonorm, 1)
macros["QcritFloorAn"] = f(Q_crit_floor, 1)
macros["pfloor"] = f(P_FLOOR, 2)
macros["gsum"] = f(g_s, 1)
macros["Msum"] = f(Mk.sum(), 1)
for m, tag in (("cond", "CD"), ("cond_nonorm", "CN"), ("cond_gate_col", "CG")):
    macros[f"dyn{tag}bot"] = f(dyn_cond[m][-1, 2]); macros[f"dyn{tag}top"] = f(dyn_cond[m][-1, 0])

with open(os.path.join(HERE, "numbers.tex"), "w") as fh:
    for key, val in macros.items():
        fh.write(f"\\newcommand{{\\{key}}}{{{val}}}\n")

print("A", A, "L", Lk, "M", Mk, "W", W)
print("phiw_ref", phiw_ref, "phiF_ref", phiF_ref, "q2d", res2d["q"], "qnode", qf_ref, "qm", qm_ref)
print("T_leg", T_leg, "T/L", T_leg / Lk, "dp", dp_leg)
print("A: phiw", phiw_A, "q", qfA + qmA)
print("C: phiw", phiw_C, "qf", qfC, "qm", qmC, "Qcrit", Q_crit)
print("seed: s", s_seed, "q", qf_seed + qm_seed, "Qcrit", Q_crit_seed)
for m, v in rem.items():
    print("rem", m, "phiw", v[0], "qf", v[1], "qm", v[2], "T", v[3])
for m, v in dyn.items():
    print("dyn", m, v[-1], v[0])
for m in cond_modes:
    print("cond", m, "phiw", cond_fp[m][0], "q", rates(m), "alpha", cond_alpha[m], "share", cond_fracshare[m], "Qcrit", cond_Qcrit[m])
print("Qcrit nonorm an", Q_crit_nonorm, "floor an", Q_crit_floor)
for m, v in dyn_cond.items():
    print("dyncond", m, v[-1], v[0])
print("pw_qb0", pw_qb0, "wmax mm", res2d["wmax"] * 1e3, "feed", res2d["feed"])
