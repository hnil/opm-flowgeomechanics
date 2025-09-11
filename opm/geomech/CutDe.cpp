/*
This is a modified version of code from https://github.com/tbenthompson/cutde
in the file cytde/common.cu.
Only the interface code is modified to easier call it from c++
*/
/*
# Modified from code presented in:
# Nikkhoo, M., Walter, T. R. (2015): Triangular dislocation: an analytical,
# artefact-free solution. - Geophysical Journal International, 201,
# 1117-1139. doi: 10.1093/gji/ggv035

# Original documentation:
# TDdispFS
# calculates displacements associated with a triangular dislocation in an
# elastic full-space.
#
# TD: Triangular Dislocation
# EFCS: Earth-Fixed Coordinate System
# TDCS: Triangular Dislocation Coordinate System
# ADCS: Angular Dislocation Coordinate System
#
# INPUTS
# X, Y and Z:
# Coordinates of calculation points in EFCS (East, North, Up). X, Y and Z
# must have the same size.
#
# P1,P2 and P3:
# Coordinates of TD vertices in EFCS.
#
# Ss, Ds and Ts:
# TD slip vector components (Strike-slip, Dip-slip, Tensile-slip).
#
# nu:
# Poisson's ratio.
#
# OUTPUTS
# ue, un and uv:
# Calculated displacement vector components in EFCS. ue, un and uv have
# the same unit as Ss, Ds and Ts in the inputs.
#
# Original documentation license:
# Copyright (c) 2014 Mehdi Nikkhoo
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files
# (the "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to permit
# persons to whom the Software is furnished to do so, subject to the
# following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
# NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
# USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <opm/geomech/CutDe.hpp>

#include <array>
#include <cmath>

#define WITHIN_KERNEL /*;*/

#ifndef M_PI
#define M_PI (3.14159265358979323846264338327950288)
#endif

#ifndef EPS
#define EPS (1e-14)
#endif

namespace
{
ddm::Real3
add3(const ddm::Real3& a, const ddm::Real3& b)
{
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

ddm::Real6
add6(const ddm::Real6& a, const ddm::Real6& b)
{
    return {a.x + b.x, a.y + b.y, a.z + b.z, a.a + b.a, a.b + b.b, a.c + b.c};
}

ddm::Real3
sub3(const ddm::Real3& a, const ddm::Real3& b)
{
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

ddm::Real3
div_scalar3(const ddm::Real3& a, const ddm::Real b)
{
    return {a.x / b, a.y / b, a.z / b};
}

ddm::Real3
mul_scalar3(const ddm::Real3& a, const ddm::Real b)
{
    return {a.x * b, a.y * b, a.z * b};
}

ddm::Real2
make2(const ddm::Real x, const ddm::Real y)
{
    return {x, y};
}

} // Anonymous namespace

ddm::Real3
ddm::make3(const Real x, const Real y, const Real z)
{
    return {x, y, z};
}

ddm::Real6
ddm::make6(const Real x, const Real y, const Real z, const Real a, const Real b, const Real c)
{
    return {x,y,z,a,b,c};
}

namespace
{

ddm::Real
dot2(const ddm::Real2& a, const ddm::Real2& b)
{
    return a.x * b.x + a.y * b.y;
}

ddm::Real
dot3(const ddm::Real3& a, const ddm::Real3& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

ddm::Real
length3(const ddm::Real3& a)
{
    return std::hypot(a.x, a.y, a.z);
}

ddm::Real3
normalize3(const ddm::Real3& a)
{
    return div_scalar3(a, length3(a));
}

ddm::Real3
negate3(const ddm::Real3 a)
{
    return div_scalar3(a, -1.0);
}

WITHIN_KERNEL ddm::Real2
transform2(const ddm::Real2& a, const ddm::Real2 b, const ddm::Real2 v)
{
    return {dot2(a, v), dot2(b, v)};
}

WITHIN_KERNEL ddm::Real2
inv_transform2(const ddm::Real2& a, const ddm::Real2& b, const ddm::Real2& v)
{
    return {a.x * v.x + b.x * v.y, a.y * v.x + b.y * v.y};
}

WITHIN_KERNEL ddm::Real3
transform3(const ddm::Real3& a, const ddm::Real3& b, const ddm::Real3& c, const ddm::Real3& v)
{
    return {dot3(a, v), dot3(b, v), dot3(c, v)};
}

WITHIN_KERNEL ddm::Real3
inv_transform3(const ddm::Real3& a, const ddm::Real3& b, const ddm::Real3& c, const ddm::Real3& v)
{
    return {a.x * v.x + b.x * v.y + c.x * v.z,
            a.y * v.x + b.y * v.y + c.y * v.z,
            a.z * v.x + b.z * v.y + c.z * v.z};
}

WITHIN_KERNEL ddm::Real3
cross3(const ddm::Real3& u, const ddm::Real3& v)
{
    return {u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x};
}

WITHIN_KERNEL ddm::Real6
tensor_transform3(const ddm::Real3& a,
                  const ddm::Real3& b,
                  const ddm::Real3& c,
                  const ddm::Real6& tensor)
{
    const auto A = std::array {a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z};

    auto out = ddm::Real6 {};

    out.x = A[0] * A[0] * tensor.x + 2 * A[0] * A[3] * tensor.a + 2 * A[0] * A[6] * tensor.b
        + 2 * A[3] * A[6] * tensor.c + A[3] * A[3] * tensor.y + A[6] * A[6] * tensor.z;

    out.y = A[1] * A[1] * tensor.x + 2 * A[1] * A[4] * tensor.a + 2 * A[1] * A[7] * tensor.b
        + 2 * A[4] * A[7] * tensor.c + A[4] * A[4] * tensor.y + A[7] * A[7] * tensor.z;

    out.z = A[2] * A[2] * tensor.x + 2 * A[2] * A[5] * tensor.a + 2 * A[2] * A[8] * tensor.b
        + 2 * A[5] * A[8] * tensor.c + A[5] * A[5] * tensor.y + A[8] * A[8] * tensor.z;

    out.a = A[0] * A[1] * tensor.x + (A[0] * A[4] + A[1] * A[3]) * tensor.a
        + (A[0] * A[7] + A[1] * A[6]) * tensor.b + (A[7] * A[3] + A[6] * A[4]) * tensor.c
        + A[4] * A[3] * tensor.y + A[6] * A[7] * tensor.z;

    out.b = A[0] * A[2] * tensor.x + (A[0] * A[5] + A[2] * A[3]) * tensor.a
        + (A[0] * A[8] + A[2] * A[6]) * tensor.b + (A[8] * A[3] + A[6] * A[5]) * tensor.c
        + A[5] * A[3] * tensor.y + A[6] * A[8] * tensor.z;

    out.c = A[1] * A[2] * tensor.x + (A[2] * A[4] + A[1] * A[5]) * tensor.a
        + (A[2] * A[7] + A[1] * A[8]) * tensor.b + (A[7] * A[5] + A[8] * A[4]) * tensor.c
        + A[4] * A[5] * tensor.y + A[7] * A[8] * tensor.z;

    return out;
}

WITHIN_KERNEL int
trimodefinder(const ddm::Real3& obs,
              const ddm::Real3& tri0,
              const ddm::Real3& tri1,
              const ddm::Real3& tri2)
{
    // trimodefinder calculates the normalized barycentric coordinates of
    // the points with respect to the TD vertices and specifies the appropriate
    // artefact-free configuration of the angular dislocations for the
    // calculations. The input matrices x, y and z share the same size and
    // correspond to the y, z and x coordinates in the TDCS, respectively. p1,
    // p2 and p3 are two-component matrices representing the y and z coordinates
    // of the TD vertices in the TDCS, respectively.
    // The components of the output (trimode) corresponding to each calculation
    // points, are 1 for the first configuration, -1 for the second
    // configuration and 0 for the calculation point that lie on the TD sides.

    const auto a = ((tri1.z - tri2.z) * (obs.y - tri2.y) + (tri2.y - tri1.y) * (obs.z - tri2.z))
        / ((tri1.z - tri2.z) * (tri0.y - tri2.y) + (tri2.y - tri1.y) * (tri0.z - tri2.z));

    const auto b = ((tri2.z - tri0.z) * (obs.y - tri2.y) + (tri0.y - tri2.y) * (obs.z - tri2.z))
        / ((tri1.z - tri2.z) * (tri0.y - tri2.y) + (tri2.y - tri1.y) * (tri0.z - tri2.z));

    const auto c = 1 - a - b;

    int result = 1;
    if ((a <= 0 && b > c && c > a) || (b <= 0 && c > a && a > b) || (c <= 0 && a > b && b > c)) {
        result = -1;
    }

    if ((a == 0 && b >= 0 && c >= 0) || (a >= 0 && b == 0 && c >= 0) || (a >= 0 && b >= 0 && c == 0)) {
        result = 0;
    }

    if (result == 0 && obs.x != 0) {
        result = 1;
    }

    return result;
}

WITHIN_KERNEL ddm::Real3
AngDisDisp(const ddm::Real x,
           const ddm::Real y,
           const ddm::Real z,
           const ddm::Real alpha,
           const ddm::Real bx,
           const ddm::Real by,
           const ddm::Real bz,
           const ddm::Real nu)
{
    const auto cosA = std::cos(alpha);
    const auto sinA = std::sin(alpha);
    const auto eta = y * cosA - z * sinA;
    const auto zeta = y * sinA + z * cosA;
    const auto r = std::hypot(x, y, z);

    const auto ux = bx / 8 / M_PI / (1 - nu) * (x * y / r / (r - z) - x * eta / r / (r - zeta));
    const auto vx = bx / 8 / M_PI / (1 - nu)
        * (eta * sinA / (r - zeta) - y * eta / r / (r - zeta) + y * y / r / (r - z)
           + (1 - 2 * nu) * (cosA * log(r - zeta) - log(r - z)));

    const auto wx = bx / 8 / M_PI / (1 - nu)
        * (eta * cosA / (r - zeta) - y / r - eta * z / r / (r - zeta)
           - (1 - 2 * nu) * sinA * log(r - zeta));

    const auto uy = by / 8 / M_PI / (1 - nu)
        * (x * x * cosA / r / (r - zeta) - x * x / r / (r - z)
           - (1 - 2 * nu) * (cosA * log(r - zeta) - log(r - z)));

    const auto vy = by * x / 8 / M_PI / (1 - nu)
        * (y * cosA / r / (r - zeta) - sinA * cosA / (r - zeta) - y / r / (r - z));

    const auto wy
        = by * x / 8 / M_PI / (1 - nu) * (z * cosA / r / (r - zeta) - cosA * cosA / (r - zeta) + 1 / r);

    const auto uz
        = bz * sinA / 8 / M_PI / (1 - nu) * ((1 - 2 * nu) * log(r - zeta) - x * x / r / (r - zeta));

    const auto vz = bz * x * sinA / 8 / M_PI / (1 - nu) * (sinA / (r - zeta) - y / r / (r - zeta));
    const auto wz = bz * x * sinA / 8 / M_PI / (1 - nu) * (cosA / (r - zeta) - z / r / (r - zeta));

    return ddm::make3(ux + uy + uz, vx + vy + vz, wx + wy + wz);
}

WITHIN_KERNEL ddm::Real3
TDSetupD(const ddm::Real3& obs,
         const ddm::Real alpha,
         const ddm::Real3& slip,
         const ddm::Real nu,
         const ddm::Real3& TriVertex,
         const ddm::Real3& SideVec)
{
    // TDSetupD transforms coordinates of the calculation points as well as
    // slip vector components from ADCS into TDCS. It then calculates the
    // displacements in ADCS and transforms them into TDCS.

    const auto A1 = make2(SideVec.z, -SideVec.y);
    const auto A2 = make2(SideVec.y, SideVec.z);
    const auto r1 = transform2(A1, A2, make2(obs.y - TriVertex.y, obs.z - TriVertex.z));

    const auto y1 = r1.x;
    const auto z1 = r1.y;

    const auto r2 = transform2(A1, A2, make2(slip.y, slip.z));
    const auto by1 = r2.x;
    const auto bz1 = r2.y;

    const auto uvw = AngDisDisp(obs.x, y1, z1, -M_PI + alpha, slip.x, by1, bz1, nu);

    const auto r3 = inv_transform2(A1, A2, make2(uvw.y, uvw.z));
    const auto v = r3.x;
    const auto w = r3.y;

    return ddm::make3(uvw.x, v, w);
}

WITHIN_KERNEL ddm::Real6
AngDisStrain(const ddm::Real x,
             const ddm::Real y,
             const ddm::Real z,
             const ddm::Real alpha,
             const ddm::Real bx,
             const ddm::Real by,
             const ddm::Real bz,
             const ddm::Real nu)
{
    // AngDisStrain calculates the strains associated with an angular
    // dislocation in an elastic full-space.

    const auto cosA = std::cos(alpha);
    const auto sinA = std::sin(alpha);
    const auto eta = y * cosA - z * sinA;
    const auto zeta = y * sinA + z * cosA;

    const auto x2 = x * x;
    const auto y2 = y * y;
    const auto z2 = z * z;
    const auto r2 = x2 + y2 + z2;
    const auto r = sqrt(r2);
    const auto r3 = r * r2;
    const auto rz = r * (r - z);
    const auto rmz = (r - z);
    const auto r2z2 = r2 * rmz * rmz;
    const auto r3z = r3 * rmz;

    const auto W = zeta - r;
    const auto W2 = W * W;
    const auto Wr = W * r;
    const auto W2r = W2 * r;
    const auto Wr3 = W * r3;
    const auto W2r2 = W2 * r2;

    const auto C = (r * cosA - z) / Wr;
    const auto S = (r * sinA - y) / Wr;

    // Partial derivatives of the Burgers' function
    const auto rFi_rx = (eta / r / (r - zeta) - y / r / (r - z)) / 4 / M_PI;
    const auto rFi_ry = (x / r / (r - z) - cosA * x / r / (r - zeta)) / 4 / M_PI;
    const auto rFi_rz = (sinA * x / r / (r - zeta)) / 4 / M_PI;

    auto out = ddm::Real6 {};

    out.x = bx * (rFi_rx)
        + bx / 8 / M_PI / (1 - nu)
            * (eta / Wr + eta * x2 / W2r2 - eta * x2 / Wr3 + y / rz - x2 * y / r2z2 - x2 * y / r3z)
        - by * x / 8 / M_PI / (1 - nu)
            * (((2 * nu + 1) / Wr + x2 / W2r2 - x2 / Wr3) * cosA + (2 * nu + 1) / rz - x2 / r2z2
               - x2 / r3z)
        + bz * x * sinA / 8 / M_PI / (1 - nu) * ((2 * nu + 1) / Wr + x2 / W2r2 - x2 / Wr3);

    out.y = by * (rFi_ry)
        + bx / 8 / M_PI / (1 - nu)
            * ((1 / Wr + S * S - y2 / Wr3) * eta + (2 * nu + 1) * y / rz - y * y * y / r2z2
               - y * y * y / r3z - 2 * nu * cosA * S)
        - by * x / 8 / M_PI / (1 - nu)
            * (1 / rz - y2 / r2z2 - y2 / r3z + (1 / Wr + S * S - y2 / Wr3) * cosA)
        + bz * x * sinA / 8 / M_PI / (1 - nu) * (1 / Wr + S * S - y2 / Wr3);

    out.z = bz * (rFi_rz)
        + bx / 8 / M_PI / (1 - nu)
            * (eta / W / r + eta * C * C - eta * z2 / Wr3 + y * z / r3 + 2 * nu * sinA * C)
        - by * x / 8 / M_PI / (1 - nu) * ((1 / Wr + C * C - z2 / Wr3) * cosA + z / r3)
        + bz * x * sinA / 8 / M_PI / (1 - nu) * (1 / Wr + C * C - z2 / Wr3);

    out.a = bx * (rFi_ry) / 2 + by * (rFi_rx) / 2
        - bx / 8 / M_PI / (1 - nu)
            * (x * y2 / r2z2 - nu * x / rz + x * y2 / r3z - nu * x * cosA / Wr + eta * x * S / Wr
               + eta * x * y / Wr3)
        + by / 8 / M_PI / (1 - nu)
            * (x2 * y / r2z2 - nu * y / rz + x2 * y / r3z + nu * cosA * S + x2 * y * cosA / Wr3
               + x2 * cosA * S / Wr)
        - bz * sinA / 8 / M_PI / (1 - nu) * (nu * S + x2 * S / Wr + x2 * y / Wr3);

    out.b = bx * (rFi_rz) / 2 + bz * (rFi_rx) / 2
        - bx / 8 / M_PI / (1 - nu)
            * (-x * y / r3 + nu * x * sinA / Wr + eta * x * C / Wr + eta * x * z / Wr3)
        + by / 8 / M_PI / (1 - nu)
            * (-x2 / r3 + nu / r + nu * cosA * C + x2 * z * cosA / Wr3 + x2 * cosA * C / Wr)
        - bz * sinA / 8 / M_PI / (1 - nu) * (nu * C + x2 * C / Wr + x2 * z / Wr3);

    out.c = by * (rFi_rz) / 2 + bz * (rFi_ry) / 2
        + bx / 8 / M_PI / (1 - nu)
            * (y2 / r3 - nu / r - nu * cosA * C + nu * sinA * S + eta * sinA * cosA / W2
               - eta * (y * cosA + z * sinA) / W2r + eta * y * z / W2r2 - eta * y * z / Wr3)
        - by * x / 8 / M_PI / (1 - nu)
            * (y / r3 + sinA * cosA * cosA / W2 - cosA * (y * cosA + z * sinA) / W2r
               + y * z * cosA / W2r2 - y * z * cosA / Wr3)
        - bz * x * sinA / 8 / M_PI / (1 - nu)
            * (y * z / Wr3 - sinA * cosA / W2 + (y * cosA + z * sinA) / W2r - y * z / W2r2);

    return out;
}

WITHIN_KERNEL ddm::Real6
TDSetupS(const ddm::Real3& obs,
         const ddm::Real alpha,
         const ddm::Real3& slip,
         const ddm::Real nu,
         const ddm::Real3& TriVertex,
         const ddm::Real3& SideVec)
{
    // TDSetupS transforms coordinates of the calculation points as well as
    // slip vector components from ADCS into TDCS. It then calculates the
    // strains in ADCS and transforms them into TDCS.

    // Transformation matrix
    const auto A1 = make2(SideVec.z, -SideVec.y);
    const auto A2 = make2(SideVec.y, SideVec.z);

    // Transform coordinates of the calculation points from TDCS into ADCS
    const auto r1 = transform2(A1, A2, make2(obs.y - TriVertex.y, obs.z - TriVertex.z));
    const auto y1 = r1.x;
    const auto z1 = r1.y;

    // Transform the in-plane slip vector components from TDCS into ADCS
    const auto r2 = transform2(A1, A2, make2(slip.y, slip.z));
    const auto by1 = r2.x;
    const auto bz1 = r2.y;

    // Calculate strains associated with an angular dislocation in ADCS
    const auto out_adcs = AngDisStrain(obs.x, y1, z1, -M_PI + alpha, slip.x, by1, bz1, nu);

    // Transform strains from ADCS into TDCS
    const auto B0 = ddm::make3(1.0, 0.0, 0.0);
    const auto B1 = ddm::make3(0.0, A1.x, A1.y);
    const auto B2 = ddm::make3(0.0, A2.x, A2.y);

    return tensor_transform3(B0, B1, B2, out_adcs);
}

// only used for haspace
void
setup_tde(ddm::Real3& transformed_obs,
          std::array<ddm::Real3, 3>& transformed_tri,
          ddm::Real& A,
          ddm::Real& B,
          ddm::Real& C,
          ddm::Real3& e12,
          ddm::Real3& e13,
          ddm::Real3& e23,
          ddm::Real3& Vnorm,
          ddm::Real3& Vstrike,
          ddm::Real3& Vdip,
          int& mode,
          const ddm::Real3& obs,
          const std::array<ddm::Real3, 3>& tri_prefix,
          const ddm::Real3& slip,
          const ddm::Real& nu,
          const bool is_halfspace)
{
    // printf("is_halfspace: %s\n", ${is_halfspace} ? "true":"false");
    // Real3
    Vnorm = normalize3(cross3(sub3(tri_prefix[1], tri_prefix[0]), sub3(tri_prefix[2], tri_prefix[0])));

    const auto eY = ddm::make3(0.0f, 1.0f, 0.0f);
    const auto eZ = ddm::make3(0.0f, 0.0f, 1.0f);

    // Real3
    Vstrike = cross3(eZ, Vnorm);
    if (length3(Vstrike) == 0) {
        Vstrike = mul_scalar3(eY, Vnorm.z);
        if (is_halfspace) {
            // For horizontal elements in case of half-space calculation!!!
            // Correct the strike vector of image dislocation only
            if (tri_prefix[0].z > 0) {
                Vstrike = negate3(Vstrike);
            }
        }
    }

    Vstrike = normalize3(Vstrike);

    // Real3
    Vdip = cross3(Vnorm, Vstrike);

    // Real3
    transformed_obs = transform3(Vnorm, Vstrike, Vdip, sub3(obs, tri_prefix[1]));

    // Real3
    transformed_tri[0] = transform3(Vnorm, Vstrike, Vdip, sub3(tri_prefix[0], tri_prefix[1]));

    // Real3
    transformed_tri[1] = ddm::make3(0.0f, 0.0f, 0.0f);

    // Real3
    transformed_tri[2] = transform3(Vnorm, Vstrike, Vdip, sub3(tri_prefix[2], tri_prefix[1]));

    e12 = normalize3(sub3(transformed_tri[1], transformed_tri[0]));
    e13 = normalize3(sub3(transformed_tri[2], transformed_tri[0]));
    e23 = normalize3(sub3(transformed_tri[2], transformed_tri[1]));

    A = std::acos(dot3(e12, e13));
    B = std::acos(dot3(negate3(e12), e23));
    C = std::acos(dot3(e23, e13));

    mode = trimodefinder(transformed_obs, transformed_tri[0], transformed_tri[1], transformed_tri[2]);
}

} // Anonymous namespace

ddm::Real3
ddm::disp_fs(const Real3& obs, const std::array<Real3, 3>& tri, const Real3& slip, const Real nu)
{
    Real3 transformed_obs;
    std::array<Real3, 3> transformed_tri;
    Real A, B, C;
    Real3 e12, e13, e23;
    Real3 Vnorm, Vstrike, Vdip;
    int mode;

    setup_tde(transformed_obs,
              transformed_tri,
              A,
              B,
              C,
              e12,
              e13,
              e23,
              Vnorm,
              Vstrike,
              Vdip,
              mode,
              obs,
              tri,
              slip,
              nu,
              /*is_halfspace*/ false);

    Real3 out {};
    if (mode == 1) {
        // Calculate first angular dislocation contribution
        const auto r1Tp = TDSetupD(transformed_obs, A, slip, nu, transformed_tri[0], negate3(e13));

        // Calculate second angular dislocation contribution
        const auto r2Tp = TDSetupD(transformed_obs, B, slip, nu, transformed_tri[1], e12);

        // Calculate third angular dislocation contribution
        const auto r3Tp = TDSetupD(transformed_obs, C, slip, nu, transformed_tri[2], e23);

        out = add3(add3(r1Tp, r2Tp), r3Tp);
    } else if (mode == -1) {
        // Calculate first angular dislocation contribution
        const auto r1Tn = TDSetupD(transformed_obs, A, slip, nu, transformed_tri[0], e13);

        // Calculate second angular dislocation contribution
        const auto r2Tn = TDSetupD(transformed_obs, B, slip, nu, transformed_tri[1], negate3(e12));

        // Calculate third angular dislocation contribution
        const auto r3Tn = TDSetupD(transformed_obs, C, slip, nu, transformed_tri[2], negate3(e23));

        out = add3(add3(r1Tn, r2Tn), r3Tn);
    } else {
        out = make3(NAN, NAN, NAN);
    }

    const auto a = make3(-transformed_obs.x,
                         transformed_tri[0].y - transformed_obs.y,
                         transformed_tri[0].z - transformed_obs.z);

    const auto b = negate3(transformed_obs);

    const auto c = make3(-transformed_obs.x,
                         transformed_tri[2].y - transformed_obs.y,
                         transformed_tri[2].z - transformed_obs.z);

    const auto na = length3(a);
    const auto nb = length3(b);
    const auto nc = length3(c);

    const auto FiN
        = a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x) + a.z * (b.x * c.y - b.y * c.x);

    const auto FiD = na * nb * nc + dot3(a, b) * nc + dot3(a, c) * nb + dot3(b, c) * na;
    const auto Fi = -2 * std::atan2(FiN, FiD) / (4 / M_PI);

    // Calculate the complete displacement vector components in TDCS
    out = add3(out, mul_scalar3(slip, Fi));

    // Transform the complete displacement vector components from TDCS into EFCS
    return inv_transform3(Vnorm, Vstrike, Vdip, out);
}

ddm::Real6
ddm::strain_fs(const Real3& obs, const std::array<Real3, 3>& tri, const Real3& slip, const Real nu)
{
    Real3 transformed_obs;
    std::array<Real3, 3> transformed_tri;
    Real A, B, C;
    Real3 e12, e13, e23;
    Real3 Vnorm, Vstrike, Vdip;
    int mode;

    setup_tde(transformed_obs,
              transformed_tri,
              A,
              B,
              C,
              e12,
              e13,
              e23,
              Vnorm,
              Vstrike,
              Vdip,
              mode,
              obs,
              tri,
              slip,
              nu,
              /*is_halfspace*/ false);

    Real6 out;
    if (mode == 1) {
        // Calculate first angular dislocation contribution
        Real6 comp1 = TDSetupS(transformed_obs, A, slip, nu, transformed_tri[0], negate3(e13));
        // Calculate second angular dislocation contribution
        Real6 comp2 = TDSetupS(transformed_obs, B, slip, nu, transformed_tri[1], e12);
        // Calculate third angular dislocation contribution
        Real6 comp3 = TDSetupS(transformed_obs, C, slip, nu, transformed_tri[2], e23);

        out = add6(add6(comp1, comp2), comp3);
    } else if (mode == -1) {
        // Calculate first angular dislocation contribution
        const Real6 comp1 = TDSetupS(transformed_obs, A, slip, nu, transformed_tri[0], e13);

        // Calculate second angular dislocation contribution
        const Real6 comp2 = TDSetupS(transformed_obs, B, slip, nu, transformed_tri[1], negate3(e12));

        // Calculate third angular dislocation contribution
        const Real6 comp3 = TDSetupS(transformed_obs, C, slip, nu, transformed_tri[2], negate3(e23));

        out = add6(add6(comp1, comp2), comp3);
    } else {
        out = make6(NAN, NAN, NAN, NAN, NAN, NAN);
    }

    return tensor_transform3(Vnorm, Vstrike, Vdip, out);
}
