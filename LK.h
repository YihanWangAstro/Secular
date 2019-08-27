#ifndef QUAD_LK_H
#define QUAD_LK_H

#include<algorithm>
#include "tools.h"
#include "SpaceHub/src/orbits/orbits.hpp"

namespace secular {

    class BasicConst {
    public:
        BasicConst(double _m1, double _m2, double _m3)
                : m1_{_m1},
                  m2_{_m2},
                  m3_{_m3},
                  m12_{_m1 + _m2},
                  m_tot_{_m1 + _m2 + _m3},
                  mu_in_{_m1 * _m2 / (_m1 + _m2)},
                  mu_out_{(_m1 + _m2) * _m3 / (_m1 + _m2 + _m3)},
                  a_in_coef_{calc_a_coef(m12_, mu_in_)},
                  a_out_coef_{calc_a_coef(m_tot_, mu_out_)},
                  SA_acc_coef_{consts::G * _m3 / mu_out_} {}

        BasicConst() = default;

        READ_GETTER(double, m1, m1_);

        READ_GETTER(double, m2, m2_);

        READ_GETTER(double, m3, m3_);

        READ_GETTER(double, m12, m12_);

        READ_GETTER(double, m_tot, m_tot_);

        READ_GETTER(double, mu_in, mu_in_);

        READ_GETTER(double, mu_out, mu_out_);

        READ_GETTER(double, a_in_coef, a_in_coef_);

        READ_GETTER(double, a_out_coef, a_out_coef_);

        READ_GETTER(double, SA_acc_coef, SA_acc_coef_);
    private:
        double m1_{0};
        double m2_{0};
        double m3_{0};
        double m12_{0};
        double m_tot_{0};
        double mu_in_{0};
        double mu_out_{0};
        double a_in_coef_{0};
        double a_out_coef_{0};
        double SA_acc_coef_{0};

        inline double calc_a_coef(double m, double mu) {
            return 1 / (consts::G * m) / mu / mu;
        }
    };

    double t_k_quad(double m_in, double m_out, double a_in, double a_out_eff) {
        double ratio = a_out_eff / sqrt(a_in);
        return consts::r_G_sqrt * sqrt(m_in) / m_out * ratio * ratio * ratio;
    }

    double normed_oct_epsilon(double m1, double m2, double a_in, double a_out, double c_out_sqr) {
        return fabs(m1 - m2) / (m1 + m2) * a_in / a_out / c_out_sqr;
    }

    auto unpack_init(size_t b, std::vector<double> const &v) {
        return std::make_tuple(v[b], v[b + 1], v[b + 2], v[b + 3], v[b + 4], v[b + 5], v[b + 6], v[b + 7], v[b + 8], v[b + 9], v[b + 10], v[b + 11]);
    }

    template<typename Container>
    void initilize_orbit_args(bool DA, size_t spin_num, Container &c, std::vector<double> const &o, size_t off_set) {
        if (DA) {
            initilize_DA(spin_num, c, o, off_set);
        } else {
            initilize_SA(spin_num, c, o, off_set);
        }
    }

    template<typename Container>
    void initilize_DA(size_t spin_num, Container &c, std::vector<double> const &args, size_t off_set) {
        auto[m1, m2, m3, a_in, a_out, e_in, e_out, omega_in, omega_out, Omega_in, i_in, i_out] = unpack_init(off_set, args);

        double Omega_out = Omega_in - 180;

        deg_to_rad(omega_in, omega_out, Omega_in, Omega_out, i_in, i_out);

        auto[j1x, j1y, j1z] = secular::unit_j(i_in, Omega_in);

        double L1 = secular::calc_angular_mom(m1, m2, a_in) * sqrt(1 - e_in * e_in);

        c.set_L1(L1 * j1x, L1 * j1y, L1 * j1z);

        auto[j2x, j2y, j2z] = secular::unit_j(i_out, Omega_out);

        double L2 = secular::calc_angular_mom(m1 + m2, m3, a_out) * sqrt(1 - e_out * e_out);

        c.set_L2(L2 * j2x, L2 * j2y, L2 * j2z);

        auto[e1x, e1y, e1z] = secular::unit_e(i_in, omega_in, Omega_in);

        c.set_e1(e_in * e1x, e_in * e1y, e_in * e1z);

        auto[e2x, e2y, e2z] = secular::unit_e(i_out, omega_out, Omega_out);

        c.set_e2(e_out * e2x, e_out * e2y, e_out * e2z);

        std::copy_n(args.begin() + off_set + 12, spin_num * 3, c.spin_begin());
    }

    template<typename Container>
    void initilize_SA(size_t spin_num, Container &c, std::vector<double> const &args, size_t off_set) {
        auto[m1, m2, m3, a_in, a_out, e_in, e_out, omega_in, omega_out, Omega_in, i_in, i_out] = unpack_init(off_set, args);

        double Omega_out = Omega_in - 180;

        double M_nu = args[22];

        deg_to_rad(omega_in, omega_out, Omega_in, Omega_out, i_in, i_out, M_nu);

        auto[j1x, j1y, j1z] = secular::unit_j(i_in, Omega_in);

        double L1 = secular::calc_angular_mom(m1, m2, a_in) * sqrt(1 - e_in * e_in);

        c.set_L1(L1 * j1x, L1 * j1y, L1 * j1z);

        auto[e1x, e1y, e1z] = secular::unit_e(i_in, omega_in, Omega_in);

        c.set_e1(e_in * e1x, e_in * e1y, e_in * e1z);

        double E_nu = space::orbit::calc_eccentric_anomaly(M_nu, e_out);

        double cosE = cos(E_nu);

        double nu_out = space::orbit::calc_true_anomaly(E_nu, e_out);//acos( ( cosE - o.e_out)/ (1 - o.e_out*cosE) );

        double r = a_out * (1 - e_out * cosE);

        auto[px, py, pz] = secular::unit_e(i_out, omega_out + nu_out, Omega_out);

        c.set_r(r * px, r * py, r * pz);

        double v = sqrt(consts::G * (m1 + m2 + m3) / (a_out * (1 - e_out * e_out)));

        double ve = -v * sin(nu_out);

        double vv = v * (e_out + cos(nu_out));

        auto[e2x, e2y, e2z] = secular::unit_e(i_out, omega_out, Omega_out);

        auto[vx, vy, vz] = secular::unit_peri_v(i_out, omega_out, Omega_out);

        c.set_v(ve * e2x + vv * vx, ve * e2y + vv * vy, ve * e2z + vv * vz);

        std::copy_n(args.begin() + off_set + 13, spin_num * 3, c.spin_begin());
    }

    template<bool Oct, typename Args, typename Container>
    inline void double_aved_LK(Args const &args, Container const &var, Container &dvar, double t) {
        auto[e1_sqr, j1_sqr, j1, L1_norm, L_in, a_in] = calc_orbit_args(args.a_in_coef(), var.L1x(), var.L1y(), var.L1z(), var.e1x(), var.e1y(),
                                                                        var.e1z());

        auto[e2_sqr, j2_sqr, j2, L2_norm, L_out, a_out] = calc_orbit_args(args.a_out_coef(), var.L2x(), var.L2y(), var.L2z(), var.e2x(), var.e2y(),
                                                                          var.e2z());
        /*---------------------------------------------------------------------------*\
            unit vectors
        \*---------------------------------------------------------------------------*/
        double n1x = var.L1x() / L1_norm, n1y = var.L1y() / L1_norm, n1z = var.L1z() / L1_norm;

        double n2x = var.L2x() / L2_norm, n2y = var.L2y() / L2_norm, n2z = var.L2z() / L2_norm;
        /*---------------------------------------------------------------------------*\
            dot production
        \*---------------------------------------------------------------------------*/
        double dn1n2 = dot(n1x, n1y, n1z, n2x, n2y, n2z);

        double de1n2 = dot(var.e1x(), var.e1y(), var.e1z(), n2x, n2y, n2z);
        /*---------------------------------------------------------------------------*\
            cross production
        \*---------------------------------------------------------------------------*/
        auto const[cn1n2_x, cn1n2_y, cn1n2_z] = cross(n1x, n1y, n1z, n2x, n2y, n2z);

        auto const[cn1e1_x, cn1e1_y, cn1e1_z] = cross(n1x, n1y, n1z, var.e1x(), var.e1y(), var.e1z());

        auto const[ce1n2_x, ce1n2_y, ce1n2_z] = cross(var.e1x(), var.e1y(), var.e1z(), n2x, n2y, n2z);

        auto const[ce2n1_x, ce2n1_y, ce2n1_z] = cross(var.e2x(), var.e2y(), var.e2z(), n1x, n1y, n1z);

        auto const[ce1e2_x, ce1e2_y, ce1e2_z] = cross(var.e1x(), var.e1y(), var.e1z(), var.e2x(), var.e2y(), var.e2z());

        auto const[cn2e2_x, cn2e2_y, cn2e2_z] = cross(n2x, n2y, n2z, var.e2x(), var.e2y(), var.e2z());
        /*---------------------------------------------------------------------------*\
            combinations
        \*---------------------------------------------------------------------------*/
        double a_out_eff = a_out * j2;

        double quad_coef = 0.75 / t_k_quad(args.m1() + args.m2(), args.m3(), a_in, a_out_eff);

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

        dvar.set_L1(dLx, dLy, dLz);

        dvar.set_L2(-dLx, -dLy, -dLz);

        dvar.set_e1(B1 * ce1n2_x + B2 * cn1e1_x + B3 * cn1n2_x, B1 * ce1n2_y + B2 * cn1e1_y + B3 * cn1n2_y,
                    B1 * ce1n2_z + B2 * cn1e1_z + B3 * cn1n2_z);

        dvar.set_e2(C1 * ce1e2_x + C2 * ce2n1_x + C3 * cn2e2_x, C1 * ce1e2_y + C2 * ce2n1_y + C3 * cn2e2_y,
                    C1 * ce1e2_z + C2 * ce2n1_z + C3 * cn2e2_z);
    }

    template<bool Oct, typename Args, typename Container>
    inline void single_aved_LK(Args const &args, Container const &var, Container &dvar, double t) {
        auto[e1_sqr, j1_sqr, j1, L1_norm, L_in, a_in] = calc_orbit_args(args.a_in_coef(), var.L1x(), var.L1y(), var.L1z(), var.e1x(), var.e1y(),
                                                                        var.e1z());

        double r2 = norm2(var.rx(), var.ry(), var.rz());

        double r = sqrt(r2);
        /*---------------------------------------------------------------------------*\
            unit vectors
        \*---------------------------------------------------------------------------*/
        double n1x = var.L1x() / L1_norm, n1y = var.L1y() / L1_norm, n1z = var.L1z() / L1_norm;

        double rhox = var.rx() / r, rhoy = var.ry() / r, rhoz = var.rz() / r;
        /*---------------------------------------------------------------------------*\
            dot production
        \*---------------------------------------------------------------------------*/
        double dn1r = dot(n1x, n1y, n1z, rhox, rhoy, rhoz);

        double de1r = dot(var.e1x(), var.e1y(), var.e1z(), rhox, rhoy, rhoz);
        /*---------------------------------------------------------------------------*\
            cross production
        \*---------------------------------------------------------------------------*/
        auto const[cn1r_x, cn1r_y, cn1r_z] = cross(n1x, n1y, n1z, rhox, rhoy, rhoz);

        auto const[cn1e1_x, cn1e1_y, cn1e1_z] = cross(n1x, n1y, n1z, var.e1x(), var.e1y(), var.e1z());

        auto const[ce1r_x, ce1r_y, ce1r_z] = cross(var.e1x(), var.e1y(), var.e1z(), rhox, rhoy, rhoz);
        /*---------------------------------------------------------------------------*\
            combinations
        \*---------------------------------------------------------------------------*/
        double quad_coef = 1.5 / t_k_quad(args.m1() + args.m2(), args.m3(), a_in, r);

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

        double D = -0.75 * args.SA_acc_coef() * args.mu_in() * a_in * a_in;

        double acc_r = -args.SA_acc_coef() * (args.m1() + args.m2()) / r3 + D * (25 * de1r * de1r - 5 * j1_sqr * dn1r * dn1r + 1 - 6 * e1_sqr) / r5;

        double acc_n = D * 2 * j1_sqr * dn1r / r4;

        double acc_e = -D * 10 * de1r / r4;

        dvar.set_L1(A1 * ce1r_x + A2 * cn1r_x, A1 * ce1r_y + A2 * cn1r_y, A1 * ce1r_z + A2 * cn1r_z);

        dvar.set_e1(B1 * cn1r_x + B2 * ce1r_x + B3 * cn1e1_x, B1 * cn1r_y + B2 * ce1r_y + B3 * cn1e1_y, B1 * cn1r_z + B2 * ce1r_z + B3 * cn1e1_z);

        dvar.set_r(var.vx(), var.vy(), var.vz());

        dvar.set_v(acc_r * var.rx() + acc_n * n1x + acc_e * var.e1x(), acc_r * var.ry() + acc_n * n1y + acc_e * var.e1y(),
                   acc_r * var.rz() + acc_n * n1z + acc_e * var.e1z());
    }

} // namespace secular
#endif
