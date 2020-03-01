#ifndef SEECUAR_STELLAR_H
#define SECULAR_STELLAR_H

#include <cmath>
#include <tuple>

#include "SpaceHub/src/orbits/orbits.hpp"
#include "SpaceHub/src/rand-generator.hpp"
#include "tools.h"
namespace secular {

double stellar_age(double m, double Z) { return 10e10 * pow(m, -2.5); }

auto kick(double _1D_sigma) {
  double vx = space::random::Normal(0, _1D_sigma);
  double vy = space::random::Normal(0, _1D_sigma);
  double vz = space::random::Normal(0, _1D_sigma);

  return std::make_tuple(vx, vy, vz);
}

class StellarConst {
 public:
  StellarConst(double m1, double m2, double m3, double Z1 = 0, double Z2 = 0, double Z3 = 0) {
    m1_ = m1;
    m2_ = m2;
    m3_ = m3;
    Z1_ = Z1;
    Z2_ = Z2;
    Z3_ = Z3;
    m1_age_ = stellar_age(m1, Z1);
    m2_age_ = stellar_age(m2, Z2);
    m3_age_ = stellar_age(m3, Z3);
  }

  StellarConst() = default;

  READ_GETTER(double, m1_age, m1_age_);

  READ_GETTER(double, m2_age, m2_age_);

  READ_GETTER(double, m3_age, m3_age_);

  READ_GETTER(double, m1, m1_);

  READ_GETTER(double, m2, m2_);

  READ_GETTER(double, m3, m3_);

  READ_GETTER(bool, m1_dead, m1_dead_);

  READ_GETTER(bool, m2_dead, m2_dead_);

  READ_GETTER(bool, m3_dead, m3_dead_);

  void make_m1_exploded() { evolve(m1_, Z1_, m1_dead_); }

  void make_m2_exploded() { evolve(m2_, Z2_, m2_dead_); }

  void make_m3_exploded() { evolve(m3_, Z3_, m3_dead_); }

 private:
  double m1_{1};
  double m2_{1};
  double m3_{1};
  double Z1_{0};
  double Z2_{0};
  double Z3_{0};
  double m1_age_{0};
  double m2_age_{0};
  double m3_age_{0};
  bool m1_dead_{false};
  bool m2_dead_{false};
  bool m3_dead_{false};

  void evolve(double &m, double Z, bool &is_dead) {
    m = m * 0.9;
    is_dead = true;
  }
};

template <typename Args, typename Container>
auto to_pos_vel_in(Args const &args, Container const &var, double M_nu) {
  auto [e1_sqr, j1_sqr, j1, L1_norm, L_in, a_in] =
      calc_orbit_args(args.a_in_coef(), var.L1x(), var.L1y(), var.L1z(), var.e1x(), var.e1y(), var.e1z());

  double e_in = sqrt(e1_sqr);

  double E_nu = space::orbit::M_anomaly_to_E_anomaly(M_nu, e_in);

  double cosE = cos(E_nu);

  double nu_in = space::orbit::E_anomaly_to_T_anomaly(E_nu, e_in);

  double r = a_in * (1 - e_in * cosE);

  double Ltx = var.L1x() + var.L2x();

  double Lty = var.L1y() + var.L2y();

  double Ltz = var.L1z() + var.L2z();

  auto [Asx, Asy, Asz] = cross(Ltx, Lty, Ltz, var.L1x(), var.L1y(), var.L1z());

  double omega_in = angle(Asx, Asy, Asz, var.e1x(), var.e1y(), var.e1z());

  double Omega_in = 0;

  double i_in = angle(var.L1x(), var.L1y(), var.L1z(), Ltx, Lty, Ltz);

  auto [npx, npy, npz] = secular::unit_e(i_in, omega_in + nu_in, Omega_in);

  double v = sqrt(consts::G * (args.m1() + args.m2() + args.m3()) / (a_in * (1 - e_in * e_in)));

  double ve = -v * sin(nu_in);

  double vv = v * (e_in + cos(nu_in));

  auto [e2x, e2y, e2z] = secular::unit_e(i_in, omega_in, Omega_in);

  auto [nvx, nvy, nvz] = secular::unit_peri_v(i_in, omega_in, Omega_in);

  return std::make_tuple(r * npx, r * npy, r * npz, ve * e2x + vv * nvx, ve * e2y + vv * nvy, ve * e2z + vv * nvz);
}
/*
template <typename Args, typename Container>
auto to_pos_vel_out(Args const &args, Container const &var, double M_nu) {
  auto [e2_sqr, j2_sqr, j2, L2_norm, L_out, a_out] =
      calc_orbit_args(args.a_out_coef(), var.L2x(), var.L2y(), var.L2z(), var.e1x(), var.e1y(), var.e1z());

  double e_in = sqrt(e1_sqr);

  double E_nu = space::orbit::M_anomaly_to_E_anomaly(M_nu, e_in);

  double cosE = cos(E_nu);

  double nu_in = space::orbit::E_anomaly_to_T_anomaly(E_nu, e_in);

  double r = a_in * (1 - e_in * cosE);

  double Ltx = var.L1x() + var.L2x();

  double Lty = var.L1y() + var.L2y();

  double Ltz = var.L1z() + var.L2z();

  auto [Asx, Asy, Asy] = cross(Ltx, Lty, Ltz, var.L1x(), var.L1y(), var.L1z());

  double omega_in = angle(Asx, Asy, Asz, var.e1x(), var.e1y(), var.e1z());

  double Omega_in = 0;

  double i_in = angle(var.L1x(), var.L1y(), var.L1z(), Ltx, Lty, Ltz);

  auto [npx, npy, npz] = secular::unit_e(i_in, omega_in + nu_in, Omega_in);

  double v = sqrt(consts::G * (args.m1() + args.m2() + args.m3()) / (a_in * (1 - e_in * e_in)));

  double ve = -v * sin(nu_in);

  double vv = v * (e_in + cos(nu_in));

  auto [e2x, e2y, e2z] = secular::unit_e(i_in, omega_in, Omega_in);

  auto [nvx, nvy, nvz] = secular::unit_peri_v(i_in, omega_in, Omega_in);

  return std::make_tuple(r * npx, r * npy, r * npz, ve * e2x + vv * nvx, ve * e2y + vv * nvy, ve * e2z + vv * nvz)
}*/

template <typename Ctrl, typename Args, typename Container>
bool stellar(Ctrl const &ctrl, Args const &args, Container const &var, Container &dvar, double t) {
  if (!args.m1_dead() && t > args.m1_age()) {
    double M_anomaly = space::random::Uniform(-consts::pi, consts::pi);

    auto [px, py, pz, vx, vy, vz] = to_pos_vel(args, var, M_anomaly);

    args.make_m1_exploded();

    auto [kick_vx, kick_vy, kick_vz] = kick(1);

    vx += kick_vx, vy += kick_vy, vz += kick_vz;

    auto [ex, ey, ez] =
        space::orbit::calc_runge_lenz_vector(consts::G * (args.m1() + args.m2()), px, py, pz, vx, vy, vz);

    if (norm(ex, ey, ez) >= 1) {
      return false;
    }

    auto [Lx, Ly, Lz] = cross_with_coef(args.mu_in(), px, py, pz, vx, vy, vz);

    dvar.add_e1(ex - var.e1x(), ey - var.e1y(), ez - var.e1z());
    dvar.add_L1(Lx - var.L1x(), Ly - var.L1y(), Lz - var.L1z());
  }

  return true;
}
}  // namespace secular
#endif