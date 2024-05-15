#pragma once
namespace ddm{
typedef double Real;
typedef struct Real2 {
    Real x;
    Real y;
} Real2;

typedef struct Real3 {
    Real x;
    Real y;
    Real z;
} Real3;

typedef struct Real6 {
    Real x;
    Real y;
    Real z;
    Real a;
    Real b;
    Real c;
} Real6;

Real3 make3(double x, double y, double z);
Real6 make6(double x, double y, double z, double a, double b, double c);
Real3 disp_fs(Real3 obs, std::array<Real3,3> tri, Real3 slip, Real nu);
Real6 strain_fs(Real3 obs, std::array<Real3,3> tri, Real3 slip, Real nu);

}
