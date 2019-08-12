#ifndef SECULAR_TOOL_H
#define SECULAR_TOOL_H

#include <cmath>
#include "vector3d.h"

namespace secular{
namespace consts{
    constexpr double pi = 3.14159265358979323;
    constexpr double G = 4 * pi * pi;
    constexpr double C = 6.32397263e4;
    constexpr double r_G_sqrt = 1.0 / (2 * pi);

    constexpr double year = 1;
    }

double calc_angular_mom(double m_in, double m_out, double a)
{
    double mu = m_in * m_out / (m_in + m_out);
    return mu * sqrt(consts::G * (m_in + m_out) * a);
}

Vec3d unit_j(double i, double Omega)
{
    double sini = sin(i);
    double cosi = cos(i);
    double sinO = sin(Omega);
    double cosO = cos(Omega);
    return Vec3d{sini * sinO, -sini * cosO, cosi};
}

Vec3d unit_e(double i, double omega, double Omega)
{
    double sini = sin(i);
    double cosi = cos(i);
    double sinO = sin(Omega);
    double cosO = cos(Omega);
    double sino = sin(omega);
    double coso = cos(omega);

    return Vec3d{coso * cosO - sino * cosi * sinO, coso * sinO + sino * cosi * cosO, sino * sini};
}

Vec3d unit_peri_v(double i, double omega, double Omega)
{
    return cross(unit_j(i,Omega), unit_e(i, omega, Omega));
} 

template <typename T, size_t len>
std::ostream &operator<<(std::ostream &os, std::array<T, len> const &arr)
{
    for (auto const &c : arr)
    {
        os << c << ' ';
    }
    return os;
}

template <typename Obt>
void deg_to_rad(Obt &args)
{
    constexpr double rad = consts::pi / 180.0;
    args.omega_in *= rad;
    args.omega_out *= rad;
    args.Omega_in *= rad;
    args.Omega_out *= rad;
    args.i_in *= rad;
    args.i_out *= rad;
    args.M_nu *= rad;
}

inline double norm2(double x, double y, double z){
    return x*x + y*y + z*z;
}

inline double norm(double x, double y, double z){
    return sqrt(norm2(x,y,z));
}

inline double dot(double x1, double y1, double z1, double x2, double y2, double z2) {
    return x1*x2 + y1*y2 + z1*z2;
}

inline auto cross(double x1, double y1, double z1, double x2, double y2, double z2){
    return std::make_tuple(y1*z2 - y2*z1, z1*x2 -z2*x1, x1*y2-x2*y1);
}

inline auto calc_orbit_args(double Coef, double lx, double ly, double lz, double ex, double ey, double ez) {
    double e_sqr = norm2(ex, ey, ez);

    double j_sqr = 1 - e_sqr;

    double j = sqrt(j_sqr);

    double L_norm = norm(lx, ly, lz);

    double L = L_norm / j;

    double a = Coef * L * L;

    return std::make_tuple(e_sqr, j_sqr, j, L_norm, L, a);
}
}
#endif
