#ifndef RELATIVISTIC_H
#define RELATIVISTIC_H

#include "tools.h"

namespace secular {

  class GRConst{
  public:
      GRConst(double M, double Mu) {
          constexpr double C5 =  consts::C *consts::C* consts::C* consts::C*consts::C;
          constexpr double G3 =  consts::G* consts::G*consts::G;
          GR_coef_ = 3 * consts::G /(consts::C * consts::C) * M / Mu;//3 * pow(consts::G*M, 1.5) / (consts::C * consts::C);
          GW_L_coef_ = -6.4 * pow(consts::G, 3.5) * Mu * Mu * pow(M, 2.5) / C5;
          GW_e_coef_ = -304.0 / 15 * G3 * Mu * M * M / C5 ;
      }

      GRConst() = default;

      READ_GETTER(double, GR_coef, GR_coef_);
      READ_GETTER(double, GW_L_coef, GW_L_coef_);
      READ_GETTER(double, GW_e_coef, GW_L_coef_);
  private:
      double GR_coef_{0};
      double GW_L_coef_{0};
      double GW_e_coef_{0};
  };

    template<typename Args, typename Container>
    inline void GR_precession(Args const &args, Container const &var, Container &dvar, double t) {
        double a_eff = calc_a_eff(args.a_in_coef(), var.L1x(), var.L1y(), var.L1z(), var.e1x(), var.e1y(), var.e1z());

        double Omega = args.GR_coef() / (a_eff * a_eff * a_eff);

        dvar.add_e1(cross_with_coef(Omega, var.L1x(), var.L1y(), var.L1z(), var.e1x(), var.e1y(), var.e1z()));
    }

    template<typename Args, typename Container>
    inline void GW_radiation(Args const &args, Container const &var, Container &dvar, double t) {
        auto[e1_sqr, j1_sqr, j1, L1_norm, L_in, a_in] = calc_orbit_args(args.a_in_coef(), var.L1x(), var.L1y(), var.L1z(), var.e1x(), var.e1y(), var.e1z());

        double n1x = var.L1x() / L1_norm, n1y = var.L1y() / L1_norm, n1z = var.L1z() / L1_norm;

        double r_a = 1.0 / a_in;

        double r_a2 = r_a * r_a;

        double r_a4 = r_a2 * r_a2;

        double r_a7 = r_a4 * r_a2 * r_a;

        double r_j2 = 1.0 / j1_sqr;

        double r_j4 = r_j2 * r_j2;

        double r_j5 = r_j4 / j1;

        double GW_L_coef = args.GW_L_coef() * sqrt(r_a7) * r_j4 * (1 + 0.875 * e1_sqr);

        double GW_e_coef = args.GW_e_coef() * r_a4 * r_j5 * (1 + 121.0 / 304 * e1_sqr);

        dvar.add_L1(GW_L_coef * n1x, GW_L_coef * n1y, GW_L_coef * n1z);

        dvar.add_e1(GW_e_coef * e1x, GW_e_coef * e1y, GW_e_coef * e1z);
    }
} // namespace secular
#endif
