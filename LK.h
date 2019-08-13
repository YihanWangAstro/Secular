#ifndef QUAD_LK_H
#define QUAD_LK_H

#include "tools.h"
#include "SpaceHub/src/orbits/orbits.hpp"

namespace secular {

    double t_k_quad(double m_in, double m_out, double a_in, double a_out_eff) {
        double ratio = a_out_eff / sqrt(a_in);
        return consts::r_G_sqrt * sqrt(m_in) / m_out * ratio * ratio * ratio;
    }

    double normed_oct_epsilon(double m1, double m2, double a_in, double a_out, double c_out_sqr) {
        return fabs(m1 - m2) / (m1 + m2) * a_in / a_out / c_out_sqr;
    }

    template<typename Container, typename OrbitArg>
    void initilize_DA(Container &c, OrbitArg const &o) {
        Vec3d L1 = secular::calc_angular_mom(o.m1, o.m2, o.a_in) * sqrt(1 - o.e_in * o.e_in) *
                   secular::unit_j(o.i_in, o.Omega_in);

        Vec3d L2 = secular::calc_angular_mom(o.m1 + o.m2, o.m3, o.a_out) * sqrt(1 - o.e_out * o.e_out) *
                   secular::unit_j(o.i_out, o.Omega_out);

        Vec3d e1 = o.e_in * secular::unit_e(o.i_in, o.omega_in, o.Omega_in);

        Vec3d e2 = o.e_out * secular::unit_e(o.i_out, o.omega_out, o.Omega_out);

        c[0] = L1.x, c[1] = L1.y, c[2] = L1.z;

        c[3] = e1.x, c[4] = e1.y, c[5] = e1.z;

        c[6] = L2.x, c[7] = L2.y, c[8] = L2.z;

        c[9] = e2.x, c[10] = e2.y, c[11] = e2.z;
    }

    template<typename Container, typename OrbitArg>
    void initilize_SA(Container &c, OrbitArg const &o) {
        Vec3d L1 = secular::calc_angular_mom(o.m1, o.m2, o.a_in) * sqrt(1 - o.e_in * o.e_in) *
                   secular::unit_j(o.i_in, o.Omega_in);

        Vec3d e1 = o.e_in * secular::unit_e(o.i_in, o.omega_in, o.Omega_in);

        double E_nu = space::orbit::calc_eccentric_anomaly(o.M_nu, o.e_out);

        double cosE = cos(E_nu);

        double nu_out = space::orbit::calc_true_anomaly(E_nu, o.e_out);//acos( ( cosE - o.e_out)/ (1 - o.e_out*cosE) );

        double r = o.a_out * (1 - o.e_out * cosE);

        Vec3d r_out = r * secular::unit_e(o.i_out, o.omega_out + nu_out, o.Omega_out);

        double v = sqrt(consts::G * (o.m1 + o.m2 + o.m3) / (o.a_out * (1 - o.e_out * o.e_out)));

        Vec3d v_out = v * (sin(nu_out) * secular::unit_e(o.i_out, o.omega_out, o.Omega_out) +
                           (o.e_out + cos(nu_out)) * secular::unit_peri_v(o.i_out, o.omega_out, o.Omega_out));

        c[0] = L1.x, c[1] = L1.y, c[2] = L1.z;

        c[3] = e1.x, c[4] = e1.y, c[5] = e1.z;

        c[6] = r_out.x, c[7] = r_out.y, c[8] = r_out.z;

        c[9] = v_out.x, c[10] = v_out.y, c[11] = v_out.z;

        std::cout << r_out << ' ' << v_out << "\n";
    }

    template<bool Oct, typename Args, typename Container>
    inline void double_aved_LK(Args const &args, Container const &var, Container &ddt, double t) {
        /*---------------------------------------------------------------------------*\
            mapping alias
        \*---------------------------------------------------------------------------*/
        const auto[L1x, L1y, L1z] = std::tie(var[0], var[1], var[2]);

        const auto[e1x, e1y, e1z] = std::tie(var[3], var[4], var[5]);

        const auto[L2x, L2y, L2z] = std::tie(var[6], var[7], var[8]);

        const auto[e2x, e2y, e2z] = std::tie(var[9], var[10], var[11]);

        auto[dL1x, dL1y, dL1z] = std::tie(ddt[0], ddt[1], ddt[2]);

        auto[de1x, de1y, de1z] = std::tie(ddt[3], ddt[4], ddt[5]);

        auto[dL2x, dL2y, dL2z] = std::tie(ddt[6], ddt[7], ddt[8]);

        auto[de2x, de2y, de2z] = std::tie(ddt[9], ddt[10], ddt[11]);
        /*---------------------------------------------------------------------------*\
            orbital parameters calculation
        \*---------------------------------------------------------------------------*/
        auto[e1_sqr, j1_sqr, j1, L1_norm, L_in, a_in] = calc_orbit_args(args.a_in_coef, L1x, L1y, L1z, e1x, e1y, e1z);

        auto[e2_sqr, j2_sqr, j2, L2_norm, L_out, a_out] = calc_orbit_args(args.a_out_coef, L2x, L2y, L2z, e2x, e2y, e2z);
        /*---------------------------------------------------------------------------*\
            unit vectors
        \*---------------------------------------------------------------------------*/
        double n1x = L1x / L1_norm, n1y = L1y / L1_norm, n1z = L1z / L1_norm;

        double n2x = L2x / L2_norm, n2y = L2y / L2_norm, n2z = L2z / L2_norm;
        /*---------------------------------------------------------------------------*\
            dot production
        \*---------------------------------------------------------------------------*/
        double dn1n2 = dot(n1x, n1y, n1z, n2x, n2y, n2z);

        double de1n2 = dot(e1x, e1y, e1z, n2x, n2y, n2z);
        /*---------------------------------------------------------------------------*\
            cross production
        \*---------------------------------------------------------------------------*/
        auto const[cn1n2_x, cn1n2_y, cn1n2_z] = cross(n1x, n1y, n1z, n2x, n2y, n2z);

        auto const[cn1e1_x, cn1e1_y, cn1e1_z] = cross(n1x, n1y, n1z, e1x, e1y, e1z);

        auto const[ce1n2_x, ce1n2_y, ce1n2_z] = cross(e1x, e1y, e1z, n2x, n2y, n2z);

        auto const[ce2n1_x, ce2n1_y, ce2n1_z] = cross(e2x, e2y, e2z, n1x, n1y, n1z);

        auto const[ce1e2_x, ce1e2_y, ce1e2_z] = cross(e1x, e1y, e1z, e2x, e2y, e2z);

        auto const[cn2e2_x, cn2e2_y, cn2e2_z] = cross(n2x, n2y, n2z, e2x, e2y, e2z);
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

        double dLx = A1 * cn1n2_x + A2 * ce1n2_x;

        double dLy = A1 * cn1n2_y + A2 * ce1n2_y;

        double dLz = A1 * cn1n2_z + A2 * ce1n2_z;

        dL1x = dLx, dL1y = dLy, dL1z = dLz;

        dL2x = -dLx, dL2y = -dLy, dL2z = -dLz;

        de1x = B1 * ce1n2_x + B2 * cn1e1_x + B3 * cn1n2_x;

        de1y = B1 * ce1n2_y + B2 * cn1e1_y + B3 * cn1n2_y;

        de1z = B1 * ce1n2_z + B2 * cn1e1_z + B3 * cn1n2_z;

        de2x = C1 * ce1e2_x + C2 * ce2n1_x + C3 * cn2e2_x;

        de2y = C1 * ce1e2_y + C2 * ce2n1_y + C3 * cn2e2_y;

        de2z = C1 * ce1e2_z + C2 * ce2n1_z + C3 * cn2e2_z;
    }

    template<bool Oct, typename Args, typename Container>
    inline void single_aved_LK(Args const &args, Container const &var, Container &ddt, double t) {
        /*---------------------------------------------------------------------------*\
            mapping alias
        \*---------------------------------------------------------------------------*/
        const auto[L1x, L1y, L1z] = std::tie(var[0], var[1], var[2]);

        const auto[e1x, e1y, e1z] = std::tie(var[3], var[4], var[5]);

        const auto[rx, ry, rz] = std::tie(var[6], var[7], var[8]);

        const auto[vx, vy, vz] = std::tie(var[9], var[10], var[11]);

        auto[dL1x, dL1y, dL1z] = std::tie(ddt[0], ddt[1], ddt[2]);

        auto[de1x, de1y, de1z] = std::tie(ddt[3], ddt[4], ddt[5]);

        auto[drx, dry, drz] = std::tie(ddt[6], ddt[7], ddt[8]);

        auto[dvx, dvy, dvz] = std::tie(ddt[9], ddt[10], ddt[11]);
        /*---------------------------------------------------------------------------*\
            orbital parameters calculation
        \*---------------------------------------------------------------------------*/
        auto[e1_sqr, j1_sqr, j1, L1_norm, L_in, a_in] = calc_orbit_args(args.a_in_coef, L1x, L1y, L1z, e1x, e1y, e1z);

        double r2 = norm2(rx, ry, rz);

        double r = sqrt(r2);
        /*---------------------------------------------------------------------------*\
            unit vectors
        \*---------------------------------------------------------------------------*/
        double n1x = L1x / L1_norm, n1y = L1y / L1_norm, n1z = L1z / L1_norm;

        double rhox = rx / r, rhoy = ry / r, rhoz = rz / r;
        /*---------------------------------------------------------------------------*\
            dot production
        \*---------------------------------------------------------------------------*/
        double dn1r = dot(n1x, n1y, n1z, rhox, rhoy, rhoz);

        double de1r = dot(e1x, e1y, e1z, rhox, rhoy, rhoz);
        /*---------------------------------------------------------------------------*\
            cross production
        \*---------------------------------------------------------------------------*/
        auto const[cn1r_x, cn1r_y, cn1r_z] = cross(n1x, n1y, n1z, rhox, rhoy, rhoz);

        auto const[cn1e1_x, cn1e1_y, cn1e1_z] = cross(n1x, n1y, n1z, e1x, e1y, e1z);

        auto const[ce1r_x, ce1r_y, ce1r_z] = cross(e1x, e1y, e1z, rhox, rhoy, rhoz);
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

        double r3 = r2 * r;

        double r4 = r2 * r2;

        double r5 = r2 * r3;

        double D = -0.75 * args.SA_acc_coef * args.mu1 * a_in * a_in;

        double acc_r = -args.SA_acc_coef * (args.m1 + args.m2) / r3 +
                       D * (25 * de1r * de1r - 5 * j1_sqr * dn1r * dn1r + 1 - 6 * e1_sqr) / r5;

        double acc_n = D * 2 * j1_sqr * dn1r / r4;

        double acc_e = -D * 10 * de1r / r4;

        dL1x = A1 * ce1r_x + A2 * cn1r_x;

        dL1y = A1 * ce1r_y + A2 * cn1r_y;

        dL1z = A1 * ce1r_z + A2 * cn1r_z;

        de1x = B1 * cn1r_x + B2 * ce1r_x + B3 * cn1e1_x;

        de1y = B1 * cn1r_y + B2 * ce1r_y + B3 * cn1e1_y;

        de1z = B1 * cn1r_z + B2 * ce1r_z + B3 * cn1e1_z;

        drx = vx, dry = vy, drz = vz;

        dvx = acc_r * rx + acc_n * n1x + acc_e * e1x;

        dvy = acc_r * ry + acc_n * n1y + acc_e * e1y;

        dvz = acc_r * rz + acc_n * n1z + acc_e * e1z;
    }

} // namespace secular
#endif
