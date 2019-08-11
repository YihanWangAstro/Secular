#ifndef SECULAR_H
#define SECULAR_H

#include <iostream>
#include <array>
#include <tuple>
#include <fstream>
#include <memory>
#include <cmath>
#include <functional>
#include <iomanip>
#include "vector3d.h"
#include <boost/numeric/odeint.hpp>
#include "SpaceHub/src/multi-thread/multi-thread.hpp"
#include "SpaceHub/src/dev-tools.hpp"
#include "SpaceHub/src/orbits/orbits.hpp"

namespace secular
{

constexpr double pi = 3.14159265358979323;
constexpr double G = 4 * pi * pi;
constexpr double C = 6.32397263e4;
constexpr double r_G_sqrt = 1.0 / (2 * pi);

constexpr double year = 1;

double calc_angular_mom(double m_in, double m_out, double a)
{
    double mu = m_in * m_out / (m_in + m_out);
    return mu * sqrt(G * (m_in + m_out) * a);
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

double t_k_quad(double m_in, double m_out, double a_in, double a_out_eff)
{
    double ratio = a_out_eff/ sqrt(a_in);
    return r_G_sqrt * sqrt(m_in) / m_out * ratio * ratio * ratio;
}

double normed_oct_epsilon(double m1, double m2, double a_in, double a_out, double c_out_sqr) {
    return fabs(m1-m2)/(m1+m2)*a_in/a_out/c_out_sqr;
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

template <size_t SpinNum_>
struct OrbitArgs
{
    static constexpr size_t SpinNum{SpinNum_};
    double m1;
    double m2;
    double m3;
    double a_in;
    double a_out;
    double e_in;
    double e_out;
    double omega_in;
    double omega_out;
    double Omega_in;
    double Omega_out;
    double i_in;
    double i_out;
    double M_nu;
    std::array<Vec3d, SpinNum> s;

    friend std::istream &operator>>(std::istream &is, OrbitArgs &t)
    {
        is >> t.m1 >> t.m2 >> t.m3 >> t.a_in >> t.a_out >> t.e_in >> t.e_out >> t.omega_in >> t.omega_out >> t.Omega_in >> t.i_in >> t.i_out >> t.M_nu;
        t.Omega_out = t.Omega_in - 180.0;
        for (auto &ss : t.s)
        {
            is >> ss;
        }
        return is;
    }
    friend std::ostream &operator<<(std::ostream &os, OrbitArgs const &t)
    {
        os << t.m1 << ' ' << t.m2 << ' ' << t.m3 << ' ' << t.a_in << ' ' << t.a_out << ' ' << t.e_in << ' ' << t.e_out << ' ' << t.omega_in << ' ' << t.omega_out << ' ' << t.Omega_in << ' ' << t.Omega_out << ' ' << t.i_in << ' ' << t.i_out << ' ' << t.M_nu;
        for (auto const &ss : t.s)
        {
            os << ' ' << ss;
        }
        return os;
    }
};

template <typename Obt>
void deg_to_rad(Obt &args)
{
    constexpr double rad = pi / 180.0;
    args.omega_in *= rad;
    args.omega_out *= rad;
    args.Omega_in *= rad;
    args.Omega_out *= rad;
    args.i_in *= rad;
    args.i_out *= rad;
}

struct Controler
{
    size_t id;
    bool write_traj;
    bool write_end;
    double end_time;
    double dt_out;
    bool Oct;
    bool GR;
    bool GW;
    bool SL_out;
    bool LL_couple;

    friend std::istream &operator>>(std::istream &is, Controler &t)
    {
        double tmp{0};
        is >> tmp; t.id = static_cast<size_t>(tmp);
        is >> tmp; t.write_traj = static_cast<bool>(tmp);
        is >> tmp; t.write_end = static_cast<bool>(tmp);
        is >> t.end_time;
        is >> t.dt_out;
        is >> tmp; t.Oct = static_cast<bool>(tmp);
        is >> tmp; t.GR = static_cast<bool>(tmp);
        is >> tmp; t.GW = static_cast<bool>(tmp);
        is >> tmp; t.SL_out = static_cast<bool>(tmp);
        is >> tmp; t.LL_couple = static_cast<bool>(tmp);

        //sis >> t.id >> t.write_traj >> t.write_end >> t.end_time >> t.dt_out >> t.Oct >> t.GR >> t.GW >> t.SL_out >> t.LL_couple;
        return is;
    }

    friend std::ostream &operator<<(std::ostream &os, Controler const &t)
    {
        os << t.id << ' ' << t.write_traj << ' ' << t.write_end << ' ' << t.end_time << ' ' << t.dt_out << ' ' << t.Oct << ' ' << t.GR << ' ' << t.GW << ' ' << t.SL_out << ' ' << t.LL_couple;
        return os;
    }
};

template <size_t spin_num>
struct Task
{
    Controler ctrl;
    OrbitArgs<spin_num> obt_args;
    friend std::istream &operator>>(std::istream &is, Task &t)
    {
        is >> t.ctrl >> t.obt_args;
        return is;
    }
    friend std::ostream &operator<<(std::ostream &os, Task const &t)
    {
        os << t.ctrl << ' ' << t.obt_args;
        return os;
    }
};

template<typename Container, typename OrbitArg>
void initilize_DA(Container& c, OrbitArg const&o){
    Vec3d L1 = secular::calc_angular_mom(o.m1, o.m2, o.a_in) * sqrt(1 - o.e_in * o.e_in) * secular::unit_j(o.i_in, o.Omega_in);
    Vec3d L2 = secular::calc_angular_mom(o.m1 + o.m2, o.m3, o.a_out) * sqrt(1 - o.e_out * o.e_out) * secular::unit_j(o.i_out, o.Omega_out);
    Vec3d e1 = o.e_in * secular::unit_e(o.i_in, o.omega_in, o.Omega_in);
    Vec3d e2 = o.e_out * secular::unit_e(o.i_out, o.omega_out, o.Omega_out);
    c[0] = L1.x, c[1] = L1.y, c[2] = L1.z;
    c[3] = e1.x, c[4] = e1.y, c[5] = e1.z;
    c[6] = L2.x, c[7] = L2.y, c[8] = L2.z;
    c[9] = e2.x, c[10] = e2.y, c[11] = e2.z;
}

template<typename Container, typename OrbitArg>
void initilize_SA(Container& c, OrbitArg const&o){
    Vec3d L1 = secular::calc_angular_mom(o.m1, o.m2, o.a_in) * sqrt(1 - o.e_in * o.e_in) * secular::unit_j(o.i_in, o.Omega_in);
    Vec3d e1 = o.e_in * secular::unit_e(o.i_in, o.omega_in, o.Omega_in);
    
    double E_nu = space::orbit::calc_eccentric_anomaly(o.M_nu, o.e_out);
    double cosE = cos(E_nu);
    double nu_out = acos( ( cosE - o.e_out)/ (1 - o.e_out*cosE) );
    double r  = o.a_out * (1 - o.e_out * cosE);
    Vec3d r_out = r*secular::unit_e(o.i_out, o.omega_out + nu_out, o.Omega_out);

    double v = sqrt(G * (o.m1 + o.m2 + o.m3)/(o.a_out*(1 - o.e_out*o.e_out) ));
    Vec3d v_out = v*(sin(nu_out) * secular::unit_e(o.i_out, o.omega_out, o.Omega_out) + (o.e_out + cos(nu_out)) * secular::unit_peri_v(o.i_out, o.omega_out, o.Omega_out));

    c[0] = L1.x, c[1] = L1.y, c[2] = L1.z;
    c[3] = e1.x, c[4] = e1.y, c[5] = e1.z;
    c[6] = r_out.x, c[7] = r_out.y, c[8] = r_out.z;
    c[9] = -v_out.x, c[10] = v_out.y, c[11] = v_out.z;

    std::cout << L1 << ' ' <<  e1 << ' '<< r_out << ' ' << v_out << "!!\n";
}

enum class StopFlag
{
    shrink
};

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


struct SecularArg{
    SecularArg(double _m1, double _m2, double _m3) : m1{_m1}, m2{_m2}, m3{_m3} {
        mu1 = m1 * m2 / (m1 + m2);
        mu2 = (m1 + m2) * m3 / (m1 + m2 + m3);
        a_in_coef = 1 / (G * (m1 + m2)) / mu1 / mu1;
        a_out_coef = 1 / (G * (m1 + m2 + m3)) / mu2 / mu2; 
    } 
public:
    double m1;
    double m2;
    double m3;
    double mu1;
    double mu2;
    double a_in_coef;
    double a_out_coef;
};



template<typename Args, typename Container>
inline void double_aved_LK(Args const& args, Container const& var, Container& ddt, double t){
    /*---------------------------------------------------------------------------*\
        mapping alias
    \*---------------------------------------------------------------------------*/
    const auto [L1x, L1y, L1z] = std::tie(var[0], var[1], var[2]);

    const auto [e1x, e1y, e1z] = std::tie(var[3], var[4], var[5]);

    const auto [L2x, L2y, L2z] = std::tie(var[6], var[7], var[8]);

    const auto [e2x, e2y, e2z] = std::tie(var[9], var[10], var[11]);

    auto [dL1x, dL1y, dL1z] = std::tie(ddt[0], ddt[1], ddt[2]);

    auto [de1x, de1y, de1z] = std::tie(ddt[3], ddt[4], ddt[5]);

    auto [dL2x, dL2y, dL2z] = std::tie(ddt[6], ddt[7], ddt[8]);

    auto [de2x, de2y, de2z] = std::tie(ddt[9], ddt[10], ddt[11]);
    
    /*---------------------------------------------------------------------------*\
        orbital parameters calculation
    \*---------------------------------------------------------------------------*/
    double e1_sqr = norm2(e1x, e1y, e1z);

    double e2_sqr = norm2(e2x, e2y, e2z);

    double j1_sqr = 1 - e1_sqr;

    double j2_sqr = 1 - e2_sqr;

    double j1 = sqrt(j1_sqr);

    double j2 = sqrt(j2_sqr);

    double L1_norm = norm(L1x, L1y, L1z);

    double L2_norm = norm(L2x, L2y, L2z);

    double L_in = L1_norm / j1;

    double L_out = L2_norm / j2;

    double a_in = args.a_in_coef * L_in * L_in;

    double a_out = args.a_out_coef * L_out * L_out;
    /*---------------------------------------------------------------------------*\
        unit vectors
    \*---------------------------------------------------------------------------*/
    double n1x = L1x/L1_norm, n1y = L1y/L1_norm, n1z = L1z/L1_norm;

    double n2x = L2x/L2_norm, n2y = L2y/L2_norm, n2z = L2z/L2_norm;
    /*---------------------------------------------------------------------------*\
        dot production
    \*---------------------------------------------------------------------------*/
    double dn1n2 = dot(n1x, n1y, n1z, n2x, n2y, n2z);

    double de1n2 = dot(e1x, e1y, e1z, n2x, n2y, n2z);
    /*---------------------------------------------------------------------------*\
        cross production
    \*---------------------------------------------------------------------------*/
    auto const [cn1n2_x, cn1n2_y, cn1n2_z] = cross(n1x, n1y, n1z, n2x, n2y, n2z);

    auto const [cn1e1_x, cn1e1_y, cn1e1_z] = cross(n1x, n1y, n1z, e1x, e1y, e1z);

    auto const [ce1n2_x, ce1n2_y, ce1n2_z] = cross(e1x, e1y, e1z, n2x, n2y, n2z);

    auto const [ce2n1_x, ce2n1_y, ce2n1_z] = cross(e2x, e2y, e2z, n1x, n1y, n1z);

    auto const [ce1e2_x, ce1e2_y, ce1e2_z] = cross(e1x, e1y, e1z, e2x, e2y, e2z);

    auto const [cn2e2_x, cn2e2_y, cn2e2_z] = cross(n2x, n2y, n2z, e2x, e2y, e2z);
    /*---------------------------------------------------------------------------*\
        combinations
    \*---------------------------------------------------------------------------*/
    double a_out_eff = a_out * j2;

    double quad_coef = 0.75 / t_k_quad(args.m1 + args.m2, args.m3, a_in, a_out_eff);

    double A = quad_coef * L_in;

    double A1 = A * j1_sqr * dn1n2;

    double A2 = -A * 5 * de1n2;

    double B = quad_coef * j1;

    double B1 = B * dn1n2;

    double B2 = B * 2;

    double B3 = B * de1n2 * (-5);

    double C = quad_coef * L_in / L2_norm;

    double C1 = C * 5 * de1n2;

    double C2 = C * j1_sqr * dn1n2;

    double C3 = -C * (0.5 - 3 * e1_sqr + 12.5 * de1n2 * de1n2 - 2.5 * j1_sqr * dn1n2 * dn1n2);

    double dLx = A1*cn1n2_x + A2*ce1n2_x;

    double dLy = A1*cn1n2_y + A2*ce1n2_y;

    double dLz = A1*cn1n2_z + A2*ce1n2_z;

    dL1x = dLx; dL1y = dLy; dL1z = dLz;

    dL2x = -dLx; dL2y = -dLy; dL2z = -dLz;

    de1x = B1*ce1n2_x + B2*cn1e1_x + B3*cn1n2_x;

    de1y = B1*ce1n2_y + B2*cn1e1_y + B3*cn1n2_y;

    de1z = B1*ce1n2_z + B2*cn1e1_z + B3*cn1n2_z;

    de2x = C1*ce1e2_x + C2*ce2n1_x + C3*cn2e2_x;

    de2y = C1*ce1e2_y + C2*ce2n1_y + C3*cn2e2_y;

    de2z = C1*ce1e2_z + C2*ce2n1_z + C3*cn2e2_z;
}


template<typename Args, typename Container>
inline void single_aved_LK(Args const& args, Container const& var, Container& ddt, double t){
    /*---------------------------------------------------------------------------*\
        mapping alias
    \*---------------------------------------------------------------------------*/
    const auto [L1x, L1y, L1z] = std::tie(var[0], var[1], var[2]);

    const auto [e1x, e1y, e1z] = std::tie(var[3], var[4], var[5]);

    const auto [rx, ry, rz] = std::tie(var[6], var[7], var[8]);

    const auto [vx, vy, vz] = std::tie(var[9], var[10], var[11]);

    auto [dL1x, dL1y, dL1z] = std::tie(ddt[0], ddt[1], ddt[2]);

    auto [de1x, de1y, de1z] = std::tie(ddt[3], ddt[4], ddt[5]);

    auto [drx, dry, drz] = std::tie(ddt[6], ddt[7], ddt[8]);

    auto [dvx, dvy, dvz] = std::tie(ddt[9], ddt[10], ddt[11]);
    
    /*---------------------------------------------------------------------------*\
        orbital parameters calculation
    \*---------------------------------------------------------------------------*/
    double e1_sqr = norm2(e1x, e1y, e1z);

    double r2 = norm2(rx, ry, rz);

    double j1_sqr = 1 - e1_sqr;

    double j1 = sqrt(j1_sqr);

    double r = sqrt(r2);

    double L1_norm = norm(L1x, L1y, L1z);

    double L_in = L1_norm / j1;

    double a_in = args.a_in_coef * L_in * L_in;
    /*---------------------------------------------------------------------------*\
        unit vectors
    \*---------------------------------------------------------------------------*/
    double n1x = L1x/L1_norm, n1y = L1y/L1_norm, n1z = L1z/L1_norm;

    double rhox = rx/r, rhoy = ry/r, rhoz = rz/r;
    /*---------------------------------------------------------------------------*\
        dot production
    \*---------------------------------------------------------------------------*/
    double dn1r = dot(n1x, n1y, n1z, rhox, rhoy, rhoz);

    double de1r = dot(e1x, e1y, e1z, rhox, rhoy, rhoz);
    /*---------------------------------------------------------------------------*\
        cross production
    \*---------------------------------------------------------------------------*/
    auto const [cn1r_x, cn1r_y, cn1r_z] = cross(n1x, n1y, n1z, rhox, rhoy, rhoz);

    auto const [cn1e1_x, cn1e1_y, cn1e1_z] = cross(n1x, n1y, n1z, e1x, e1y, e1z);

    auto const [ce1r_x, ce1r_y, ce1r_z] = cross(e1x, e1y, e1z, rhox, rhoy, rhoz);
    /*---------------------------------------------------------------------------*\
        combinations
    \*---------------------------------------------------------------------------*/
    double quad_coef = 1.5 / t_k_quad(args.m1 + args.m2, args.m3, a_in, r);

    double A = quad_coef * L_in;

    double A1 = A * 5 * de1r;

    double A2 = -A * j1_sqr * dn1r;
    
    double B = quad_coef * j1;

    double B1 = B * 5 * de1r;

    double B2 = -B * dn1r;

    double B3 = -2 * B;

    double r3 = r2*r;

    double r4 = r2*r2;

    double r5 = r2*r3;

    double D = 0.75 * G * args.m3 * args.mu1 / args.mu2  * a_in * a_in;

    double acc_r = -G*args.m3/args.mu2/r3 * (args.m1 + args.m2)  - 5 * C * (5*de1r*de1r - j1_sqr*dn1r*dn1r)/r5 ;

    double acc_n = - D * 2 * j1_sqr * dn1r / r4;

    double acc_e = D * 10 * de1r / r4;

    dL1x = A1*ce1r_x + A2*cn1r_x;

    dL1y = A1*ce1r_y + A2*cn1r_y;

    dL1z = A1*ce1r_z + A2*cn1r_z;

    de1x = B1*cn1r_x + B2*ce1r_x + B3*cn1e1_x;

    de1y = B1*cn1r_y + B2*ce1r_y + B3*cn1e1_y;

    de1z = B1*cn1r_z + B2*ce1r_z + B3*cn1e1_z;

    drx = vx, dry = vy, drz = vz;

    dvx = acc_r*rx + acc_n * n1x + acc_e*e1x;

    dvy = acc_r*ry + acc_n * n1y + acc_e*e1y;

    dvz = acc_r*rz + acc_n * n1z + acc_e*e1z;
}

template<typename Args, typename Container>
inline void GR_precession(Args const& args, Container const& var, Container& ddt, double t){
    
        /*---------------------------------------------------------------------------*\
            mapping alias
        \*---------------------------------------------------------------------------*/
        const auto [L1x, L1y, L1z] = std::tie(var[0], var[1], var[2]);

        const auto [e1x, e1y, e1z] = std::tie(var[3], var[4], var[5]);

        auto [de1x, de1y, de1z] = std::tie(ddt[3], ddt[4], ddt[5]);
        /*---------------------------------------------------------------------------*\
            orbital parameters calculation
        \*---------------------------------------------------------------------------*/
        double e1_sqr = norm2(e1x, e1y, e1z);

        double j1_sqr = 1 - e1_sqr;

        double j1 = sqrt(j1_sqr);

        double L1_norm = norm(L1x, L1y, L1z);

        double L_in = L1_norm / j1;

        double a_in = args.a_in_coef * L_in * L_in;
        /*---------------------------------------------------------------------------*\
            unit vectors
        \*---------------------------------------------------------------------------*/
        double n1x = L1x/L1_norm, n1y = L1y/L1_norm, n1z = L1z/L1_norm;
        /*---------------------------------------------------------------------------*\
            cross production
        \*---------------------------------------------------------------------------*/
        auto const [cn1e1_x, cn1e1_y, cn1e1_z] = cross(n1x, n1y, n1z, e1x, e1y, e1z);
        /*---------------------------------------------------------------------------*\
            combinations
        \*---------------------------------------------------------------------------*/
        double r_a = 1.0 / a_in;

        double GR_coef = args.GR_coef * sqrt(r_a * r_a * r_a * r_a * r_a)/ j1_sqr;

        de1x += GR_coef * cn1e1_x;

        de1y += GR_coef * cn1e1_y;

        de1z += GR_coef * cn1e1_z;
}


template<typename Args, typename Container>
inline void GW_radiation(Args const& args, Container const& var, Container& ddt, double t){
        /*---------------------------------------------------------------------------*\
            mapping alias
        \*---------------------------------------------------------------------------*/
        const auto [L1x, L1y, L1z] = std::tie(var[0], var[1], var[2]);

        const auto [e1x, e1y, e1z] = std::tie(var[3], var[4], var[5]);

        auto [dL1x, dL1y, dL1z] = std::tie(ddt[0], ddt[1], ddt[2]);

        auto [de1x, de1y, de1z] = std::tie(ddt[3], ddt[4], ddt[5]);
        /*---------------------------------------------------------------------------*\
            orbital parameters calculation
        \*---------------------------------------------------------------------------*/
        double e1_sqr = norm2(e1x, e1y, e1z);

        double j1_sqr = 1 - e1_sqr;

        double j1 = sqrt(j1_sqr);

        double L1_norm = norm(L1x, L1y, L1z);

        double L_in = L1_norm / j1;

        double a_in = args.a_in_coef * L_in * L_in;
        /*---------------------------------------------------------------------------*\
            unit vectors
        \*---------------------------------------------------------------------------*/
        double n1x = L1x/L1_norm, n1y = L1y/L1_norm, n1z = L1z/L1_norm;
        /*---------------------------------------------------------------------------*\
            combinations
        \*---------------------------------------------------------------------------*/
        double r_a = 1.0 / a_in;

        double r_a2 = r_a * r_a;

        double r_a4 = r_a2 * r_a2;

        double r_a7 = r_a4 * r_a2 * r_a;

        double r_j2 = 1.0 / j1_sqr;

        double r_j4 = r_j2 * r_j2;

        double r_j5 = r_j4 / j1;

        double GW_L_coef = args.GW_L_coef * sqrt(r_a7) * r_j4 * (1 + 0.875 * e1_sqr);

        double GW_e_coef = args.GW_e_coef * r_a4 * r_j5  * (1 + 121.0/304 * e1_sqr);

        dL1x += GW_L_coef * n1x;

        dL1y += GW_L_coef * n1y;

        dL1z += GW_L_coef * n1z;

        de1x += GW_e_coef * e1x;

        de1y += GW_e_coef * e1y;

        de1z += GW_e_coef * e1z;
}

inline double dsdt_coupling_coef(double C, double a, double j_sqr){
    double r_a = 1.0/a;
    double r_a2 = r_a * r_a;
    double r_a4 = r_a2 * r_a2;
    double r_a5 = r_a4 * r_a;
    return C*sqrt(r_a5)/j_sqr;
}

template<typename Args, typename Container>
inline void spin_evolve(Args const& args, Container const& var, Container& ddt, double t){
        /*---------------------------------------------------------------------------*\
            mapping alias
        \*---------------------------------------------------------------------------*/
        const auto [L1x, L1y, L1z] = std::tie(var[0], var[1], var[2]);

        const auto [e1x, e1y, e1z] = std::tie(var[3], var[4], var[5]);

        auto [dL1x, dL1y, dL1z] = std::tie(ddt[0], ddt[1], ddt[2]);

        auto [de1x, de1y, de1z] = std::tie(ddt[3], ddt[4], ddt[5]);
        /*---------------------------------------------------------------------------*\
            orbital parameters calculation
        \*---------------------------------------------------------------------------*/
        double e1_sqr = norm2(e1x, e1y, e1z);

        double j1_sqr = 1 - e1_sqr;

        double j1 = sqrt(j1_sqr);

        double L1_norm = norm(L1x, L1y, L1z);

        double L_in = L1_norm / j1;

        double a_in = args.a_in_coef * L_in * L_in;
        /*---------------------------------------------------------------------------*\
            unit vectors
        \*---------------------------------------------------------------------------*/
        double n1x = L1x/L1_norm, n1y = L1y/L1_norm, n1z = L1z/L1_norm;

        /*---------------------------------------------------------------------------*\
            combinations
        \*---------------------------------------------------------------------------*/  
}


template<typename Args, typename ...Func>
inline auto serilize(Args const& args, Func...func) {
    return [&](auto const&x, auto &dxdt, double t) {
        (func(args, x, dxdt,t),...);
    };
}


} // namespace secular
#endif