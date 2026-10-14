#!/usr/bin/env python3
"""Crossflow example with the properties of test/m51w/TEST.DATA and formation damage (WINJDAM).

Generates figs/m51w_*.pdf and numbers_m51w.tex (macros prefixed with W) for
crossflow_example.tex.  Potentials in bar, rates in reservoir m3/d, conductances in m3/(d bar).
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

g = 9.80665
bar = 1e5
CONV = bar * 86400.0
mD = 9.869233e-16

# ------------------------------------------------ deck: test/m51w/TEST.DATA
dx = 2.0                                   # DXV/DYV/DZ = 2 m
h = 2.0
kcomp = np.arange(36, 41)                  # COMPDAT 'B-3H' 26 26 36 40
zc = 1900.0 + h * (kcomp - 0.5)            # TOPS 1900
zo = 1900.0 + h * (38 - 0.5)               # WSEED K = 38
perm = 1000 * mD                           # PERMX = PERMY, layers 26-60
rw = 0.216 / 2                             # COMPDAT diameter 0.216
B_w = 1.02643 * np.exp(-3.7876e-5 * (195 - 79))  # PVTW at ~195 bar (EQUIL)
rho_w = 1026.0 / B_w                       # DENSITY water 1026
mu = 0.39831e-3                            # PVTW viscosity [Pa s]
lam = 1.0 / mu                             # water zone: krw(Sw=1) = 1, kro = 0
Q_surf = 100.0                             # WCONINJE RATE 100 sm3/d
Q = Q_surf * B_w                           # reservoir rate
E = 14e9                                   # YMODULE 14 GPa
nu = 0.35                                  # PRATIO
sig_grad = 0.167                           # STREQUIL xx/yy gradient [bar/m]
dT = 90.0 - 30.0                           # RTEMPVD 90 C, WTEMP 30 C
alpha_T = 3.0e-5                           # THERMEXR
k_cake = 10 * mD                           # WINJDAM LINRAD 10 mD
phi_cake = 0.3                             # WINJDAM porosity
conc = 10e-6                               # WINJFCNC 10 ppm (volume fraction)
# ------------------------------------------------ fracture settings (json)
d_leak = 0.25 * dx                         # reservoir.calculate_dist, legacy factor 0.25
sides = 2
feed_radius = 1.0                          # well_source_all_perfs, well_source_radius 1.0
P_FLOOR = 1e4 / bar
# ------------------------------------------------ example choices
R = 3.0                                    # fracture covers completions 37-39
p_net = 20.0                               # net pressure [bar] -> w_max ~ 1 mm
W_over_L = 100.0
gamma = 0.03                               # depth-linear mismatch for the dynamics [bar/m]
C_tank = 2.0                               # near-well storage per layer [m3/bar]
Aq = 110.0                                 # drainage of each layer to the far field
Tv = 430.0                                 # vertical coupling of adjacent near-well layers
nday = 100
N = len(zc)

r_o = 0.198 * dx
D_peace = np.log(r_o / rw)
M0 = np.full(N, 2 * np.pi * perm * h * lam / D_peace * CONV)
A_cake_m = 2 * np.pi * rw * h              # Connection::getFilterCakeArea default
K_over_kc = perm / k_cake                  # K = Kh / L
s_mis = -gamma * (zc - zo) + 0.0


class Fracture:
    """2D fracture model (cubic law, leak-off with a filter cake in series per face)."""

    def __init__(self, n=60):
        dl = 2 * R / n
        xs = -R + dl * (np.arange(n) + 0.5)
        X, Zr = np.meshgrid(xs, xs)
        r = np.hypot(X, Zr)
        self.inside = r < R
        idx = -np.ones(X.shape, int)
        idx[self.inside] = np.arange(self.inside.sum())
        self.nc = int(self.inside.sum())
        self.X, self.Z = X, Zr + zo
        self.dA = dl * dl
        self.wmax = 8 * (1 - nu**2) / (np.pi * E) * p_net * bar * R
        width = 8 * (1 - nu**2) / (np.pi * E) * p_net * bar * np.sqrt(np.clip(R**2 - r**2, 0, None))
        width = np.maximum(width, 2e-5)
        rows, cols, vals = [], [], []
        diag = np.zeros(self.nc)
        for (di, dj) in ((0, 1), (1, 0)):
            a = self.inside[:n - di, :n - dj] & self.inside[di:, dj:]
            i1 = idx[:n - di, :n - dj][a]; i2 = idx[di:, dj:][a]
            w1 = width[:n - di, :n - dj][a]; w2 = width[di:, dj:][a]
            t = lam * CONV / (12 / w1**3 + 12 / w2**3) * 2
            rows += [i1, i2]; cols += [i2, i1]; vals += [-t, -t]
            np.add.at(diag, i1, t); np.add.at(diag, i2, t)
        self.Kflow = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                                   shape=(self.nc, self.nc)) + sp.diags(diag)
        zabs = self.Z[self.inside]
        self.layer = np.clip(np.floor((zabs - (zc[0] - h / 2)) / h).astype(int), 0, N - 1)
        x_in = X[self.inside]
        self.feed = (np.abs(x_in) < feed_radius) & (zabs > zc[0] - h / 2) & (zabs < zc[-1] + h / 2)
        self.l0 = sides * lam * perm * self.dA / d_leak * CONV
        self.W = W_over_L * self.l0 * self.nc
        self.wcell = np.where(self.feed, self.W / self.feed.sum(), 0.0)

    def lcell(self, hc):
        """Fracture::updateLeakoff with a filter cake (thickness hc over all faces)."""
        res = lam * perm * self.dA / d_leak * CONV
        cake = np.where(hc > 0, lam * k_cake * self.dA / np.maximum(hc, 1e-30) * sides * CONV, np.inf)
        return sides / (1 / res + 1 / cake)

    def L(self, hc):
        lc = self.lcell(hc)
        return np.array([lc[self.layer == k].sum() for k in range(N)])

    def solve(self, phi_layers, phi_w, hc):
        lc = self.lcell(hc)
        A = (self.Kflow + sp.diags(lc + self.wcell)).tocsc()
        rhs = lc * phi_layers[self.layer] + self.wcell * phi_w
        phi = spla.spsolve(A, rhs)
        leak = lc * (phi - phi_layers[self.layer])
        q = np.array([leak[self.layer == k].sum() for k in range(N)])
        phiF = np.array([np.average(phi[self.layer == k], weights=lc[self.layer == k])
                         if np.any(self.layer == k) else phi_w for k in range(N)])
        return phi, leak, q, phiF


frac = Fracture()
present = np.array([np.any(frac.layer == k) for k in range(N)])


def flow_step(T, M, phi, s, off, phi_far, dt):
    """Implicit Flow step: tanks with storage, far-field support and vertical coupling;
    well under rate control.  Fracture part T_k (Phi_w + off_k - Phi_k), matrix part
    M_k (Phi_w + s_k - Phi_k)."""
    nn = N + 1
    Am = np.zeros((nn, nn)); b = np.zeros(nn)
    for i in range(N):
        Am[i, i] = C_tank / dt + Aq + T[i] + M[i]
        b[i] = C_tank / dt * phi[i] + Aq * phi_far[i] + T[i] * off[i] + M[i] * s[i]
        for j in (i - 1, i + 1):
            if 0 <= j < N:
                Am[i, i] += Tv; Am[i, j] -= Tv
        Am[i, N] = -(T[i] + M[i])
    Am[N, :N] = -(T + M); Am[N, N] = (T + M).sum()
    b[N] = Q - (T * off).sum() - (M * s).sum()
    x = np.linalg.solve(Am, b)
    phi, phi_w = x[:N], x[N]
    return phi, phi_w, T * (phi_w + off - phi), M * (phi_w + s - phi)


def legacy_ctf(q, dp):
    T = np.where(np.abs(dp) > 1e-14, q / np.where(np.abs(dp) > 1e-14, dp, 1.0), 0.0)
    return np.where(T > 0, T, 0.0)


def conductivity_ctf(q, dp, Lk):
    sum_q = q.sum()
    den = (Lk * np.maximum(dp, P_FLOOR)).sum()
    alpha = min(sum_q / den, 2.0) if (sum_q > 0 and den > 0) else 1.0
    return alpha * Lk


def run(mode, damage, s, phi_far=np.zeros(N), dt=1.0, nsteps=None, m_init=None, hc_init=None):
    """Sequential coupling as in the code; cake growth as in WellFilterCake (matrix,
    LINRAD) and Fracture::updateFilterCakeProps (fracture faces, scale_filtrate)."""
    phi = np.zeros(N)
    phi_w = Q / (frac.L(np.zeros(frac.nc)).sum() + M0.sum())
    hc = np.zeros(frac.nc) if hc_init is None else hc_init.copy()
    hm = np.zeros(N)
    S = np.zeros(N) if m_init is None else D_peace / m_init - D_peace
    T_old = None
    out = dict(qf=[], qm=[], m=[], Lrel=[], phiw=[], Qcrit=[], T=[], qmod=[], phi=[])
    L0 = frac.L(np.zeros(frac.nc))
    for n in range(nday if nsteps is None else nsteps):
        m = D_peace / (D_peace + S)
        M = M0 * m
        Lk = frac.L(hc)
        if mode == "reference":
            # fracture node solved with the layers (response at the current cake)
            G = frac.solve(np.zeros(N), 1.0, hc)[2]
            H = np.column_stack([frac.solve(np.eye(N)[j], 0.0, hc)[2] for j in range(N)])
            nn = N + 1
            Am = np.zeros((nn, nn)); b = np.zeros(nn)
            for i in range(N):
                Am[i, :N] = -H[i]
                Am[i, i] += C_tank / dt + Aq + M[i]
                for j in (i - 1, i + 1):
                    if 0 <= j < N:
                        Am[i, i] += Tv; Am[i, j] -= Tv
                Am[i, N] = -(G[i] + M[i])
                b[i] = C_tank / dt * phi[i] + Aq * phi_far[i]
            Am[N, :N] = H.sum(axis=0) - M; Am[N, N] = G.sum() + M.sum(); b[N] = Q
            x = np.linalg.solve(Am, b)
            phi, phi_w = x[:N], x[N]
            _, leak, q, phiF = frac.solve(phi, phi_w, hc)
            qf, qm = q, M * (phi_w - phi)
            scale = 1.0
            T = np.where(np.abs(phi_w - phi) > 1e-14, q / (phi_w - phi), 0.0)
        else:
            _, leak, q, phiF = frac.solve(phi, phi_w, hc)
            if mode == "option3":
                T = Lk; off = phiF - phi_w
            else:
                Tn = legacy_ctf(q, phi_w - phi) if mode == "legacy" else conductivity_ctf(q, phi_w - phi, Lk)
                T = Tn if T_old is None else (Tn + 2.0 * T_old) / 3.0   # wellIndicesAvrg
                T_old = T; off = s
            phi, phi_w, qf, qm = flow_step(T, M, phi, s, off, phi_far, dt)
            # scale_filtrate: Flow's fracture rate over the fracture model's leak-off
            scale = max(0.0, qf.sum()) / q.sum() if q.sum() > 1e-12 else 0.0
            _, leak, q, phiF = frac.solve(phi, phi_w, hc)
        out["qmod"].append(q); out["phi"].append(phi.copy())
        out["qf"].append(qf); out["qm"].append(qm); out["m"].append(m.copy())
        out["Lrel"].append(Lk / np.where(L0 > 0, L0, 1.0)); out["phiw"].append(phi_w); out["T"].append(T)
        s_eff = s if mode != "reference" else np.zeros(N)
        Ttot = (T if mode != "option3" else np.zeros(N)) + M
        out["Qcrit"].append((Ttot * (s - s.min())).sum())
        out["hc"] = hc.copy()
        if damage:
            # matrix cake (WellFilterCake::updateSkinFactorsAndMultipliers, LINRAD)
            tot = qf + qm
            pos = np.maximum(tot, 0.0).sum()
            xfact = max(tot.sum(), 0.0) / pos if pos > tot.sum() and pos > 1e-12 else 1.0
            share = np.where(tot > 1e-12, np.clip(qm / np.where(tot > 1e-12, tot, 1.0), 0, 1), 0.0)
            rate = np.maximum(tot, 0.0) * conc * xfact * share
            dh = rate * dt / (A_cake_m * (1 - phi_cake))
            S += K_over_kc * np.log((rw + hm + dh) / (rw + hm))
            hm += dh
            # fracture-face cake (Fracture::updateFilterCakeProps)
            hc += np.maximum(leak, 0.0) * scale * conc * dt / (frac.dA * (1 - phi_cake))
    hc_last = out.pop("hc", hc)
    res = {k: np.array(v) for k, v in out.items()}
    res["hc"] = hc_last
    return res


def steady_min(mode, gm, m_init=None, hc_init=None):
    """Smallest connection rate in the quasi-steady state (layers equilibrated)."""
    d = run(mode, False, -gm * (zc - zo) + 0.0, dt=1e6, nsteps=40, m_init=m_init, hc_init=hc_init)
    return (d["qf"][-1] + d["qm"][-1]).min()


def gamma_crit_ss(mode, m_init=None, hc_init=None, lo=0.0, hi=0.2):
    for _ in range(14):
        mid = 0.5 * (lo + hi)
        if steady_min(mode, mid, m_init, hc_init) < 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


# =============================================================== static quantities
L0 = frac.L(np.zeros(frac.nc))
G0 = frac.solve(np.zeros(N), 1.0, np.zeros(frac.nc))[2]
phiw_ref = Q / (G0.sum() + M0.sum())
_, _, q_ref, phiF_ref = frac.solve(np.zeros(N), phiw_ref, np.zeros(frac.nc))
T_leg = legacy_ctf(q_ref, np.full(N, phiw_ref))
T_cond = conductivity_ctf(q_ref, np.full(N, phiw_ref), L0)
hyd_cell = rho_w * g * h / bar
dsig_T = E * alpha_T * dT / (1 - nu) / bar
sig_h = sig_grad * zo
p_res = 195.0 - rho_w * g * (2000 - zo) / bar
gammas = np.linspace(0, 0.03, 121)
span = zc - zc.min()


def qcrit(T, gm):
    s = -gm * (zc - zo)
    return ((T + M0) * (s - s.min())).sum()


Qc_leg = np.array([qcrit(T_leg, gm) for gm in gammas])
Qc_cond = np.array([qcrit(T_cond, gm) for gm in gammas])
gamma_crit_leg = gammas[np.argmax(Qc_leg >= Q)]
gamma_crit_cond = gammas[np.argmax(Qc_cond >= Q)]
# downward drainage to the free bottom boundary (BCPROP 6 FREE): potential gradient
area_model = (51 * dx) ** 2
beta_drain = (Q / 86400 / area_model) * mu / (0.1 * perm) / bar

# =============================================================== dynamics
modes = ("reference", "legacy", "conductivity", "option3")
runs = {(m, dmg): run(m, dmg, s_mis) for m in modes for dmg in (False, True)}
# quasi-steady critical mismatch, undamaged and with the damage after nday days
cmodes = ("legacy", "conductivity", "option3")
ref_dmg = runs[("reference", True)]
gcss = {m: gamma_crit_ss(m) for m in cmodes}
gcss_d = {m: gamma_crit_ss(m, m_init=ref_dmg["m"][-1], hc_init=ref_dmg["hc"]) for m in cmodes}
skew = {m: run(m, False, -0.01 * (zc - zo) + 0.0, dt=1e6, nsteps=40) for m in cmodes}

# =============================================================== figures
cols = plt.cm.viridis(np.linspace(0, 0.9, N))
lbl = [f"K={k}" for k in kcomp]
t = np.arange(1, nday + 1)

fig, ax = plt.subplots(figsize=(5.4, 3.6))
for k in range(N):
    ax.add_patch(Rectangle((-5, zc[k] - h / 2), 10, h, fc="#efe0bd" if k % 2 else "#f6ecd6", ec="k", lw=0.4))
    ax.text(-4.8, zc[k], f"K={kcomp[k]}, z={zc[k]:.0f}", va="center", fontsize=7)
ax.add_patch(Circle((0, zo), R, fc="tab:orange", alpha=0.35, ec="tab:orange", lw=1.2))
ax.add_patch(Rectangle((-feed_radius, zc[0] - h / 2), 2 * feed_radius, N * h, fc="none", ec="tab:red",
                       ls="--", lw=0.8))
ax.plot([0, 0], [zc[0] - h / 2 - 1, zc[-1] + h / 2 + 1], color="k", lw=2.5)
for k in range(N):
    ax.plot([-0.3, 0.3], [zc[k]] * 2, color="tab:red", lw=5, solid_capstyle="butt")
ax.text(3.2, zo - 2.6, f"R = {R:.0f} m", fontsize=7)
ax.text(1.1, zc[-1] + 1.6, "feed region (well_source_radius)", fontsize=6.5, color="tab:red")
ax.set_xlim(-5, 7); ax.set_ylim(zc[-1] + h / 2 + 1, zc[0] - h / 2 - 1)
ax.set_aspect("equal"); ax.grid(False)
ax.set_xlabel("distance in fracture plane [m]"); ax.set_ylabel("depth [m]")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "m51w_geometry.pdf")); plt.close(fig)

fig, ax = plt.subplots(figsize=(4.6, 3.0))
ax.plot(gammas * 100, Qc_leg, "k", label="legacy")
ax.plot(gammas * 100, Qc_cond, "k--", label="conductivity (defaults)")
ax.axhline(Q, color="tab:red", ls=":", label="deck rate")
ax.set_xlabel(r"mismatch gradient $\gamma$ [bar/100 m]"); ax.set_ylabel(r"$Q_{crit}$ [m$^3$/d]")
ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "m51w_qcrit.pdf")); plt.close(fig)

ref_d = runs[("reference", True)]
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.8))
for k in range(N):
    axs[0].plot(t, ref_d["m"][:, k], color=cols[k], label=lbl[k])
    if present[k]:
        axs[1].plot(t, ref_d["Lrel"][:, k], color=cols[k], label=lbl[k])
axs[0].set_title("matrix CTF multiplier (WINJDAM, LINRAD)", fontsize=8)
axs[1].set_title("fracture leak-off conductance / undamaged", fontsize=8)
for ax in axs:
    ax.set_xlabel("time [d]"); ax.set_ylim(0, 1.05)
axs[0].legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "m51w_damage.pdf")); plt.close(fig)

fig, axs = plt.subplots(2, 4, figsize=(7.6, 4.6), sharex=True, sharey=True)
titles = {"reference": "reference", "legacy": "legacy", "conductivity": "conductivity",
          "option3": "option 3"}
for c, m in enumerate(modes):
    for r, dmg in enumerate((False, True)):
        ax = axs[r, c]
        d = runs[(m, dmg)]
        for k in range(N):
            ax.plot(t, d["qf"][:, k] + d["qm"][:, k], color=cols[k], label=lbl[k])
        ax.axhline(0, color="k", lw=0.6)
        ax.set_title(f"{titles[m]}, {'WINJDAM' if dmg else 'no damage'}", fontsize=7.5)
        if r == 1:
            ax.set_xlabel("time [d]")
axs[0, 0].set_ylabel(r"connection rate [m$^3$/d]"); axs[1, 0].set_ylabel(r"connection rate [m$^3$/d]")
axs[0, 0].legend(fontsize=6)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "m51w_dynamics.pdf")); plt.close(fig)

fig, ax = plt.subplots(figsize=(4.6, 3.0))
for m, ls in (("legacy", "-"), ("conductivity", "--")):
    for dmg, c in ((False, "k"), (True, "tab:red")):
        ax.plot(t, runs[(m, dmg)]["Qcrit"], color=c, ls=ls,
                label=f"{m}, {'WINJDAM' if dmg else 'no damage'}")
ax.axhline(Q, color="gray", ls=":")
ax.set_xlabel("time [d]"); ax.set_ylabel(r"$Q_{crit}$ [m$^3$/d]"); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "m51w_qcrit_time.pdf")); plt.close(fig)

fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.8), sharey=True)
for ax, dmg in zip(axs, (False, True)):
    d = runs[("legacy", dmg)]
    xb = np.arange(N)
    ax.bar(xb - 0.2, d["qf"][-1], 0.4, color="tab:red", label="Flow: fracture part of connection")
    ax.bar(xb + 0.2, d["qmod"][-1], 0.4, color="tab:orange", label="fracture model leak-off")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(xb); ax.set_xticklabels(lbl, fontsize=7)
    ax.set_title(f"legacy, day {nday}, {'WINJDAM' if dmg else 'no damage'}", fontsize=8)
axs[0].set_ylabel(r"rate [m$^3$/d]"); axs[0].legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "m51w_flow_vs_model.pdf")); plt.close(fig)

# =============================================================== numbers
def f(x, nd=2):
    return f"{x:.{nd}f}"


mac = {
    "Wrhow": f(rho_w, 0), "Wmu": f(mu * 1e3, 3), "WBw": f(B_w, 4), "WQ": f(Q, 1), "WQs": f(Q_surf, 0),
    "WR": f(R, 0), "Wpnet": f(p_net, 0), "Wwmax": f(frac.wmax * 1e3, 2), "Wdleak": f(d_leak, 1),
    "Wrw": f(rw, 3), "Wro": f(r_o, 3), "WD": f(D_peace, 3), "WM": f(M0[0], 0),
    "WLtot": f(L0.sum(), 0), "WLmid": f(L0[2], 0), "WLnb": f(L0[1], 0),
    "Wlpa": f(sides * lam * perm / d_leak * CONV, 1),
    "Wphiw": f(phiw_ref, 4), "Whyd": f(hyd_cell, 3), "Wratio": f(hyd_cell / phiw_ref, 0),
    "WTtot": f((T_leg + M0).sum(), 0), "WTleg": f(T_leg.sum(), 0),
    "Wgcl": f(gamma_crit_leg, 4), "Wgcc": f(gamma_crit_cond, 4),
    "Wgclpct": f(gamma_crit_leg / (rho_w * g / bar) * 100, 1),
    "Wgamma": f(gamma, 3), "Wdsig": f(dsig_T, 0), "Wsig": f(sig_h, 0), "Wpres": f(p_res, 1),
    "Wbeta": f(beta_drain, 4), "WAq": f(Aq, 0), "WTv": f(Tv, 0), "WC": f(C_tank, 0),
    "WKkc": f(K_over_kc, 0), "Wconc": "10", "Wkc": "10", "Wphic": f(phi_cake, 1),
}
for m in modes:
    for dmg in (False, True):
        d = runs[(m, dmg)]
        tag = {"reference": "Ref", "legacy": "Leg", "conductivity": "Con", "option3": "Opt"}[m] + ("D" if dmg else "N")
        tot = d["qf"] + d["qm"]
        mac[f"W{tag}botfirst"] = f(tot[0, -1]); mac[f"W{tag}botlast"] = f(tot[-1, -1])
        mac[f"W{tag}topfirst"] = f(tot[0, 0]); mac[f"W{tag}toplast"] = f(tot[-1, 0])
        mac[f"W{tag}minfirst"] = f(tot[0].min()); mac[f"W{tag}minlast"] = f(tot[-1].min())
        mac[f"W{tag}phiwlast"] = f(d["phiw"][-1], 3)
        mac[f"W{tag}Qcfirst"] = f(d["Qcrit"][0], 0); mac[f"W{tag}Qclast"] = f(d["Qcrit"][-1], 0)
for m, tag in (("legacy", "Leg"), ("conductivity", "Con"), ("option3", "Opt")):
    mac[f"Wgss{tag}"] = f(gcss[m], 3); mac[f"Wgssd{tag}"] = f(gcss_d[m], 3)
    mac[f"Wgsspct{tag}"] = f(gcss[m] / (rho_w * g / bar) * 100, 0)
    tot = skew[m]["qf"][-1] + skew[m]["qm"][-1]
    mac[f"Wskew{tag}top"] = f(tot[0], 1); mac[f"Wskew{tag}bot"] = f(tot[-1], 1)
dl = runs[("legacy", True)]
mac["WLegDqfTop"] = f(dl["qf"][-1][1], 1); mac["WLegDqmodTop"] = f(dl["qmod"][-1][1], 1)
mac["WLegDTTop"] = f(dl["T"][-1][1], 0); mac["WLegDLrelTop"] = f(dl["Lrel"][-1][1], 2)
mac["WLegDLrelBot"] = f(dl["Lrel"][-1][3], 2); mac["WLegDmBot"] = f(dl["m"][-1][-1], 2)
mac["WLegDphiTop"] = f(dl["phi"][-1][1], 3); mac["WLegDphiw"] = f(dl["phiw"][-1], 3)
mac["WLzero"] = f(L0[1], 0)
dl0 = runs[("legacy", False)]
mac["WLegNqfTop"] = f(dl0["qf"][-1][1], 1); mac["WLegNqmodTop"] = f(dl0["qmod"][-1][1], 1)
mac["WLegNqfLow"] = f(dl0["qf"][-1][3], 1); mac["WLegNqmodLow"] = f(dl0["qmod"][-1][3], 1)
mac["WLegNphiTop"] = f(dl0["phi"][-1][1], 3); mac["WLegNphiw"] = f(dl0["phiw"][-1], 3)
mac["WLegDqfLow"] = f(dl["qf"][-1][3], 1); mac["WLegDqmodLow"] = f(dl["qmod"][-1][3], 1)
mac["WmTen"] = f(ref_d["m"][9, 2], 2); mac["WmHundred"] = f(ref_d["m"][-1, 2], 2)
mac["WLHundred"] = f(ref_d["Lrel"][-1, 2], 2)
with open(os.path.join(HERE, "numbers_m51w.tex"), "w") as fh:
    for k_, v in mac.items():
        fh.write(f"\\newcommand{{\\{k_}}}{{{v}}}\n")

print("legacy dmg: Flow qf", dl["qf"][-1].round(2), "model q", dl["qmod"][-1].round(2), "phi", dl["phi"][-1].round(3), "phiw", dl["phiw"][-1].round(3))
dl0 = runs[("legacy", False)]
print("legacy nodmg: Flow qf", dl0["qf"][-1].round(2), "model q", dl0["qmod"][-1].round(2), "phi", dl0["phi"][-1].round(3), "phiw", dl0["phiw"][-1].round(3))
print("gcss", gcss, "gcss_d", gcss_d)
for m, d in skew.items():
    print("skew 0.01", m, (d["qf"][-1] + d["qm"][-1]).round(2))
print("zc", zc, "zo", zo, "B_w", B_w, "rho", rho_w, "Q", Q)
print("L0", L0, "M0", M0, "W", frac.W, "wmax", frac.wmax)
print("phiw_ref", phiw_ref, "q_ref", q_ref, "phiF", phiF_ref, "T_leg", T_leg, "T_cond", T_cond)
print("hyd per cell", hyd_cell, "gamma_crit", gamma_crit_leg, gamma_crit_cond, "beta", beta_drain,
      "dsig_T", dsig_T, "sig_h", sig_h, "p_res", p_res)
for (m, dmg), d in runs.items():
    tot = d["qf"] + d["qm"]
    print(m, dmg, "first", tot[0].round(2), "last", tot[-1].round(2), "phiw", round(d["phiw"][-1], 4),
          "m_last", d["m"][-1].round(2), "Lrel", d["Lrel"][-1].round(2), "Qc", round(d["Qcrit"][0], 1), round(d["Qcrit"][-1], 1))
