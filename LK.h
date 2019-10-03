#ifndef QUAD_LK_H
#define QUAD_LK_H

#include <algorithm>
#include "SpaceHub/src/orbits/orbits.hpp"
#include "tools.h"

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

  inline double calc_a_coef(double m, double mu) { return 1 / (consts::G * m) / mu / mu; }
};

enum class LK_method { DA, SA };

size_t to_index(LK_method x) {
  if (x == LK_method::SA)
    return 0;
  else if (x == LK_method::DA)
    return 1;
  else
    return 0;
}

LK_method str_to_LK_enum(std::string const &key) {
  if (case_insens_equals(key, "DA")) {
    return LK_method::DA;
  } else if (case_insens_equals(key, "SA")) {
    return LK_method::SA;
  } else {
    throw ReturnFlag::input_err;
  }
}

double t_k_quad(double m_in, double m_out, double a_in, double a_out_eff) {
  double ratio = a_out_eff / sqrt(a_in);
  return consts::r_G_sqrt * sqrt(m_in) / m_out * ratio * ratio * ratio;
}

double normed_oct_epsilon(double m1, double m2, double a_in, double a_out_eff) {
  return fabs(m1 - m2) / (m1 + m2) * a_in / a_out_eff;
}

double normed_oct_epsilon(double m1, double m2, double a_in, double a_out, double j2_sqr) {
  return fabs(m1 - m2) / (m1 + m2) * a_in / (a_out * j2_sqr);
}

/*
    auto unpack_init(size_t b, std::vector<double> const &v) {
        return std::make_tuple(v[b], v[b + 1], v[b + 2], v[b + 3], v[b + 4], v[b + 5], v[b + 6], v[b + 7], v[b + 8], v[b
   + 9], v[b + 10], v[b + 11]);
    }*/

template <typename Container, typename Iter>
void initialize_orbit_args(LK_method method, Container &c, Iter iter) {
  auto [m1, m2, m3, a_in, a_out, e_in, e_out, omega_in, omega_out, Omega_in, i_in, i_out] = unpack_args<12>(iter);

  double Omega_out = Omega_in - 180;

  double M_nu = *(iter + 12);

  deg_to_rad(omega_in, omega_out, Omega_in, Omega_out, i_in, i_out, M_nu);

  auto [j1x, j1y, j1z] = secular::unit_j(i_in, Omega_in);

  double L1 = secular::calc_angular_mom(m1, m2, a_in) * sqrt(1 - e_in * e_in);

  c.set_L1(L1 * j1x, L1 * j1y, L1 * j1z);

  auto [e1x, e1y, e1z] = secular::unit_e(i_in, omega_in, Omega_in);

  c.set_e1(e_in * e1x, e_in * e1y, e_in * e1z);

  if (method == LK_method::DA) {
    auto [j2x, j2y, j2z] = secular::unit_j(i_out, Omega_out);

    double L2 = secular::calc_angular_mom(m1 + m2, m3, a_out) * sqrt(1 - e_out * e_out);

    c.set_L2(L2 * j2x, L2 * j2y, L2 * j2z);

    auto [e2x, e2y, e2z] = secular::unit_e(i_out, omega_out, Omega_out);

    c.set_e2(e_out * e2x, e_out * e2y, e_out * e2z);
  } else if (method == LK_method::SA) {
    double E_nu = space::orbit::calc_eccentric_anomaly(M_nu, e_out);

    double cosE = cos(E_nu);

    double nu_out = space::orbit::calc_true_anomaly(E_nu, e_out);  // acos( ( cosE - o.e_out)/ (1 - o.e_out*cosE) );

    double r = a_out * (1 - e_out * cosE);

    auto [px, py, pz] = secular::unit_e(i_out, omega_out + nu_out, Omega_out);

    c.set_r(r * px, r * py, r * pz);

    double v = sqrt(consts::G * (m1 + m2 + m3) / (a_out * (1 - e_out * e_out)));

    double ve = -v * sin(nu_out);

    double vv = v * (e_out + cos(nu_out));

    auto [e2x, e2y, e2z] = secular::unit_e(i_out, omega_out, Omega_out);

    auto [vx, vy, vz] = secular::unit_peri_v(i_out, omega_out, Omega_out);

    c.set_v(ve * e2x + vv * vx, ve * e2y + vv * vy, ve * e2z + vv * vz);
  } else {
    throw ReturnFlag::input_err;
  }
  // std::copy_n(iter + 13, 9, c.spin_begin());
}

template <typename Ctrl, typename Args, typename Container>
void double_aved_LK(Ctrl const &ctrl, Args const &args, Container const &var, Container &dvar) {
  auto [e1_sqr, j1_sqr, j1, L1_norm, L_in, a_in] =
      calc_orbit_args(args.a_in_coef(), var.L1x(), var.L1y(), var.L1z(), var.e1x(), var.e1y(), var.e1z());

  auto [e2_sqr, j2_sqr, j2, L2_norm, L_out, a_out] =
      calc_orbit_args(args.a_out_coef(), var.L2x(), var.L2y(), var.L2z(), var.e2x(), var.e2y(), var.e2z());
  /*---------------------------------------------------------------------------*\
          unit vectors
  \*---------------------------------------------------------------------------*/
  double j1x = var.L1x() / L_in, j1y = var.L1y() / L_in, j1z = var.L1z() / L_in;

  double n2x = var.L2x() / L2_norm, n2y = var.L2y() / L2_norm, n2z = var.L2z() / L2_norm;
  /*---------------------------------------------------------------------------*\
          dot production
  \*---------------------------------------------------------------------------*/
  double dj1n2 = dot(j1x, j1y, j1z, n2x, n2y, n2z);

  double de1n2 = dot(var.e1x(), var.e1y(), var.e1z(), n2x, n2y, n2z);
  /*---------------------------------------------------------------------------*\
          cross production
  \*---------------------------------------------------------------------------*/
  auto const [cj1n2_x, cj1n2_y, cj1n2_z] = cross(j1x, j1y, j1z, n2x, n2y, n2z);

  auto const [cj1e1_x, cj1e1_y, cj1e1_z] = cross(j1x, j1y, j1z, var.e1x(), var.e1y(), var.e1z());

  auto const [ce1n2_x, ce1n2_y, ce1n2_z] = cross(var.e1x(), var.e1y(), var.e1z(), n2x, n2y, n2z);

  auto const [ce2j1_x, ce2j1_y, ce2j1_z] = cross(var.e2x(), var.e2y(), var.e2z(), j1x, j1y, j1z);

  auto const [ce1e2_x, ce1e2_y, ce1e2_z] = cross(var.e1x(), var.e1y(), var.e1z(), var.e2x(), var.e2y(), var.e2z());

  auto const [cn2e2_x, cn2e2_y, cn2e2_z] = cross(n2x, n2y, n2z, var.e2x(), var.e2y(), var.e2z());
  /*---------------------------------------------------------------------------*\
          combinations
  \*---------------------------------------------------------------------------*/
  double const a_out_eff = a_out * j2;

  double const quad_coef = 0.75 / t_k_quad(args.m12(), args.m3(), a_in, a_out_eff);

  double const A = quad_coef * L_in;

  double const A1 = A * dj1n2;

  double const A2 = -A * 5 * de1n2;

  double const B = quad_coef;

  double const B1 = B * dj1n2;

  double const B2 = B * 2;

  double const B3 = B * de1n2 * (-5);

  double const C = quad_coef * L_in / L2_norm;

  double const C1 = C * 5 * de1n2;

  double const C2 = C * dj1n2;

  double const C3 = -C * (0.5 - 3 * e1_sqr + 12.5 * de1n2 * de1n2 - 2.5 * dj1n2 * dj1n2);

  double const dLx = A1 * cj1n2_x + A2 * ce1n2_x;

  double const dLy = A1 * cj1n2_y + A2 * ce1n2_y;

  double const dLz = A1 * cj1n2_z + A2 * ce1n2_z;

  dvar.add_L1(dLx, dLy, dLz);

  dvar.add_e1(B1 * ce1n2_x + B2 * cj1e1_x + B3 * cj1n2_x, B1 * ce1n2_y + B2 * cj1e1_y + B3 * cj1n2_y,
              B1 * ce1n2_z + B2 * cj1e1_z + B3 * cj1n2_z);

  dvar.sub_L2(dLx, dLy, dLz);

  dvar.add_e2(C1 * ce1e2_x + C2 * ce2j1_x + C3 * cn2e2_x, C1 * ce1e2_y + C2 * ce2j1_y + C3 * cn2e2_y,
              C1 * ce1e2_z + C2 * ce2j1_z + C3 * cn2e2_z);

  if (ctrl.Oct == true) {
    double const oct_coef = -25.0 / 16 * quad_coef * normed_oct_epsilon(args.m1(), args.m2(), a_in, a_out_eff) / j2;

    auto const [cj1e2_x, cj1e2_y, cj1e2_z] = cross(j1x, j1y, j1z, var.e2x(), var.e2y(), var.e2z());

    double const de1e2 = dot(var.e1x(), var.e1y(), var.e1z(), var.e2x(), var.e2y(), var.e2z());

    double const dj1e2 = dot(j1x, j1y, j1z, var.e2x(), var.e2y(), var.e2z());

    double const shared_C1 = 1.6 * e1_sqr - 0.2 - 7 * de1n2 * de1n2 + dj1n2 * dj1n2;

    double const E1 = oct_coef * 2 * (de1e2 * dj1n2 + de1n2 * dj1e2);

    double const E2 = oct_coef * 2 * (dj1e2 * dj1n2 - 7 * de1e2 * de1n2);

    double const E3 = oct_coef * 2 * de1n2 * dj1n2;

    double const E4 = oct_coef * shared_C1;

    double const E5 = oct_coef * 3.2 * de1e2;

    double const G = -L_in / L2_norm;

    double const G1 = G * E1;

    double const G2 = G * E2;

    double const G3 = G * j2_sqr * E3;

    double const G4 = G * j2_sqr * E4;

    double const G5 =
        -oct_coef * G * ((0.4 - 3.2 * e1_sqr) * de1e2 + 14 * de1n2 * dj1e2 * dj1n2 + 7 * de1e2 * shared_C1);

    double const oct_dLx = L_in * (E1 * cj1n2_x + E2 * ce1n2_x + E3 * cj1e2_x + E4 * ce1e2_x);

    double const oct_dLy = L_in * (E1 * cj1n2_y + E2 * ce1n2_y + E3 * cj1e2_y + E4 * ce1e2_y);

    double const oct_dLz = L_in * (E1 * cj1n2_z + E2 * ce1n2_z + E3 * cj1e2_z + E4 * ce1e2_z);

    dvar.add_L1(oct_dLx, oct_dLy, oct_dLz);

    dvar.sub_L2(oct_dLx, oct_dLx, oct_dLz);

    dvar.add_e1(E1 * ce1n2_x + E2 * cj1n2_x + E3 * ce1e2_x + E4 * cj1e2_x + E5 * cj1e1_x,
                E1 * ce1n2_y + E2 * cj1n2_y + E3 * ce1e2_y + E4 * cj1e2_y + E5 * cj1e1_y,
                E1 * ce1n2_z + E2 * cj1n2_z + E3 * ce1e2_z + E4 * cj1e2_z + E5 * cj1e1_z);

    dvar.add_e2(G1 * cj1e2_x + G2 * ce1e2_x + G3 * cj1n2_x + G4 * ce1n2_x + G5 * cn2e2_x,
                G1 * cj1e2_y + G2 * ce1e2_y + G3 * cj1n2_y + G4 * ce1n2_y + G5 * cn2e2_y,
                G1 * cj1e2_z + G2 * ce1e2_z + G3 * cj1n2_z + G4 * ce1n2_z + G5 * cn2e2_z);
  }
}

template <typename Ctrl, typename Args, typename Container>
void single_aved_LK(Ctrl const &ctrl, Args const &args, Container const &var, Container &dvar) {
  auto [e1_sqr, j1_sqr, j1, L1_norm, L_in, a_in] =
      calc_orbit_args(args.a_in_coef(), var.L1x(), var.L1y(), var.L1z(), var.e1x(), var.e1y(), var.e1z());

  double const r2 = norm2(var.rx(), var.ry(), var.rz());

  double const r = sqrt(r2);
  /*---------------------------------------------------------------------------*\
          unit vectors
      \*---------------------------------------------------------------------------*/
  double const j1x = var.L1x() / L_in, j1y = var.L1y() / L_in, j1z = var.L1z() / L_in;

  double const rhox = var.rx() / r, rhoy = var.ry() / r, rhoz = var.rz() / r;
  /*---------------------------------------------------------------------------*\
          dot production
      \*---------------------------------------------------------------------------*/
  double const dj1rho = dot(j1x, j1y, j1z, rhox, rhoy, rhoz);

  double const de1rho = dot(var.e1x(), var.e1y(), var.e1z(), rhox, rhoy, rhoz);
  /*---------------------------------------------------------------------------*\
          cross production
      \*---------------------------------------------------------------------------*/
  auto const [cj1rho_x, cj1rho_y, cj1rho_z] = cross(j1x, j1y, j1z, rhox, rhoy, rhoz);

  auto const [cj1e1_x, cj1e1_y, cj1e1_z] = cross(j1x, j1y, j1z, var.e1x(), var.e1y(), var.e1z());

  auto const [ce1rho_x, ce1rho_y, ce1rho_z] = cross(var.e1x(), var.e1y(), var.e1z(), rhox, rhoy, rhoz);
  /*---------------------------------------------------------------------------*\
          combinations
      \*---------------------------------------------------------------------------*/
  double const quad_coef = 1.5 / t_k_quad(args.m12(), args.m3(), a_in, r);

  double const B1 = quad_coef * 5 * de1rho;

  double const B2 = -quad_coef * dj1rho;

  double const B3 = -2 * quad_coef;

  double const A1 = B1 * L_in;

  double const A2 = B2 * L_in;

  double const r3 = r2 * r;

  double const r4 = r2 * r2;

  double const r5 = r2 * r3;

  double const D = -0.75 * args.SA_acc_coef() * args.mu_in() * a_in * a_in;

  double const acc_r =
      -args.SA_acc_coef() * args.m12() / r3 + D * (25 * de1rho * de1rho - 5 * dj1rho * dj1rho + 1 - 6 * e1_sqr) / r5;

  double const acc_n = D * 2 * dj1rho / r4;

  double const acc_e = -D * 10 * de1rho / r4;

  dvar.add_L1(A1 * ce1rho_x + A2 * cj1rho_x, A1 * ce1rho_y + A2 * cj1rho_y, A1 * ce1rho_z + A2 * cj1rho_z);

  dvar.add_e1(B1 * cj1rho_x + B2 * ce1rho_x + B3 * cj1e1_x, B1 * cj1rho_y + B2 * ce1rho_y + B3 * cj1e1_y,
              B1 * cj1rho_z + B2 * ce1rho_z + B3 * cj1e1_z);

  dvar.set_r(var.vx(), var.vy(), var.vz());

  dvar.add_v(acc_r * var.rx() + acc_n * j1x + acc_e * var.e1x(), acc_r * var.ry() + acc_n * j1y + acc_e * var.e1y(),
             acc_r * var.rz() + acc_n * j1z + acc_e * var.e1z());

  if (ctrl.Oct == true) {
    double const epsilon = normed_oct_epsilon(args.m1(), args.m2(), a_in, r);

    double const oct_coef = quad_coef * epsilon * 5.0 / 8.0;

    double const E = 8 * e1_sqr - 1;

    double const F1 = oct_coef * 10 * dj1rho * de1rho;

    double const F2 = oct_coef * (E + 5 * dj1rho * dj1rho - 35 * de1rho * de1rho);

    double const F3 = oct_coef * 10 * dj1rho * de1rho;

    double const H1 = F1 * L_in;

    double const H2 = F2 * L_in;

    dvar.add_L1(H1 * cj1e1_x + H2 * ce1rho_x, H1 * cj1e1_y + H2 * ce1rho_y, H1 * cj1e1_z + H2 * ce1rho_z);

    dvar.add_e1(F1 * cj1e1_x + F2 * cj1rho_x + F3 * ce1rho_x, F1 * cj1e1_y + F2 * cj1rho_y + F3 * ce1rho_y,
                F1 * cj1e1_z + F2 * cj1rho_z + F3 * ce1rho_z);
  }
}

template <typename Ctrl, typename Args, typename Container>
inline void Lidov_Kozai(Ctrl const &ctrl, Args const &args, Container const &var, Container &dvar) {
  if (ctrl.ave_method == LK_method::DA) {
    double_aved_LK(ctrl, args, var, dvar);
  } else if (ctrl.ave_method == LK_method::SA) {
    single_aved_LK(ctrl, args, var, dvar);
  }
}

}  // namespace secular
#endif
