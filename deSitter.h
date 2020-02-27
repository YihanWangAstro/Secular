#ifndef DESITTER_H
#define DESITTER_H

#include "tools.h"

namespace secular {

class SLConst {
 public:
  SLConst(double m1, double m2, double m3) {
    SL_[0] = deSitter_coef(m1, m2);
    SL_[1] = deSitter_coef(m1 + m2, m3);

    SL_[2] = deSitter_coef(m2, m1);
    SL_[3] = deSitter_coef(m1 + m2, m3);

    SL_[4] = deSitter_coef(m3, m1 + m2);
    SL_[5] = SL_[4];

    LL_ = deSitter_coef(m1 + m2, m3);
  }

  SLConst() = default;

  READ_GETTER(double, LL, LL_);

  READ_GETTER(double, S1L1, SL_[0]);

  READ_GETTER(double, S1L2, SL_[1]);

  READ_GETTER(double, S2L1, SL_[2]);

  READ_GETTER(double, S2L2, SL_[3]);

  READ_GETTER(double, S3L1, SL_[4]);

  READ_GETTER(double, S3L2, SL_[5]);

 private:
  double SL_[6];
  double LL_;

  inline double deSitter_coef(double m_self, double m_other) {
    return 0.5 * consts::G / (consts::C * consts::C) * (4 + 3 * m_other / m_self);
  };
};

enum class deS { off, on, bc };

struct SLstat {
  deS Sin_Lin{deS::off};
  deS Sin_Lout{deS::off};
  deS Sout_Lin{deS::off};
  deS Sout_Lout{deS::off};
  deS LL{deS::off};
};

inline bool is_Lin_coupled(SLstat const &stat) {
  return stat.LL != deS::off || stat.Sin_Lin != deS::off || stat.Sout_Lin != deS::off;
};

inline bool is_Lout_coupled(SLstat const &stat) {
  return stat.LL != deS::off || stat.Sin_Lout != deS::off || stat.Sout_Lout != deS::off;
};

template <typename Args, typename Container>
class deSitter_arg {
 public:
  deSitter_arg(bool single_ave, SLstat const &stat, Args const &args, Container const &var) {
    const bool Lin_coupled = is_Lin_coupled(stat);
    const bool Lout_coupled = is_Lout_coupled(stat);

    if (var.spin_num > 0 && Lin_coupled) {
      a_in_eff_ = calc_a_eff(args.a_in_coef(), var.L1x(), var.L1y(), var.L1z(), var.e1x(), var.e1y(), var.e1z());

      a_in_eff3_ = a_in_eff_ * a_in_eff_ * a_in_eff_;
    }

    if (var.spin_num > 0 && Lout_coupled) {
      if (!single_ave) {
        L2x_ = var.L2x(), L2y_ = var.L2y(), L2z_ = var.L2z();

        a_out_eff_ = calc_a_eff(args.a_out_coef(), var.L2x(), var.L2y(), var.L2z(), var.e2x(), var.e2y(), var.e2z());

        a_out_eff3_ = a_out_eff_ * a_out_eff_ * a_out_eff_;
      } else {
        std::tie(L2x_, L2y_, L2z_) =
            cross_with_coef(args.mu_out(), var.rx(), var.ry(), var.rz(), var.vx(), var.vy(), var.vz());

        a_out_eff_ = norm(var.rx(), var.ry(), var.rz());

        a_out_eff3_ = a_out_eff_ * a_out_eff_ * a_out_eff_;
      }
    }

    if (var.spin_num > 0) {
      if (stat.Sin_Lin != deS::off) Omega_[0] = args.S1L1() / a_in_eff3_;
      if (stat.Sin_Lout != deS::off) Omega_[1] = args.S1L2() / a_out_eff3_;
    }

    if (var.spin_num > 1) {
      if (stat.Sin_Lin != deS::off) Omega_[2] = args.S2L1() / a_in_eff3_;
      if (stat.Sin_Lout != deS::off) Omega_[3] = args.S2L2() / a_out_eff3_;
    }

    if (var.spin_num > 2) {
      if (stat.Sout_Lin != deS::off) Omega_[4] = args.S3L1() / a_in_eff3_;
      if (stat.Sout_Lout != deS::off) Omega_[5] = args.S3L2() / a_out_eff3_;
    }

    if (stat.LL != deS::off) {
      LL_ = args.LL() / a_out_eff3_;
    }
  }

  READ_GETTER(double, LL, LL_);

  READ_GETTER(double, S1L1_Omega, Omega_[0]);

  READ_GETTER(double, S1L2_Omega, Omega_[1]);

  READ_GETTER(double, S2L1_Omega, Omega_[2]);

  READ_GETTER(double, S2L2_Omega, Omega_[3]);

  READ_GETTER(double, S3L1_Omega, Omega_[4]);

  READ_GETTER(double, S3L2_Omega, Omega_[5]);

  READ_GETTER(double, L2x, L2x_);

  READ_GETTER(double, L2y, L2y_);

  READ_GETTER(double, L2z, L2z_);

 private:
  double L2x_;
  double L2y_;
  double L2z_;
  double Omega_[3 * 2 + 1];
  double LL_;
  double a_in_eff_;
  double a_in_eff3_;
  double a_out_eff_;
  double a_out_eff3_;
};

inline auto deSitter_e_vec(double S1x, double S1y, double S1z, double Lx, double Ly, double Lz) {
  double dot_part = 3 * dot(Lx, Ly, Lz, S1x, S1y, S1z) / norm2(Lx, Ly, Lz);
  return std::make_tuple(S1x - dot_part * Lx, S1y - dot_part * Ly, S1z - dot_part * Lz);
}

template <typename Container>
auto SA_back_reaction(double Omega, double Sx, double Sy, double Sz, Container const &var) {
  auto [crvx, crvy, crvz] = cross(var.rx(), var.ry(), var.rz(), var.vx(), var.vy(), var.vz());

  double r2 = norm2(var.rx(), var.ry(), var.rz());

  double const acc_coef = Omega / r2;

  auto [csvx, csvy, csvz] = cross(Sx, Sy, Sz, var.vx(), var.vy(), var.vz());

  auto [csrx, csry, csrz] = cross(Sx, Sy, Sz, var.rx(), var.ry(), var.rz());

  double dvr = dot(var.rx(), var.ry(), var.rz(), var.vx(), var.vy(), var.vz());

  double tri_dot = dot(Sx, Sy, Sz, crvx, crvy, crvz);

  double acc_x = acc_coef * (3 * tri_dot * var.rx() + 2 * r2 * csvx - 3 * dvr * csrx);

  double acc_y = acc_coef * (3 * tri_dot * var.ry() + 2 * r2 * csvy - 3 * dvr * csry);

  double acc_z = acc_coef * (3 * tri_dot * var.rz() + 2 * r2 * csvz - 3 * dvr * csrz);

  return std::make_tuple(acc_x, acc_y, acc_z);
}

#define EVOLVE_SEL_IN(C, S)                                                                                       \
  if (C != deS::off) {                                                                                            \
    auto [dx, dy, dz] =                                                                                           \
        cross_with_coef(d.S##L1_Omega(), var.L1x(), var.L1y(), var.L1z(), var.S##x(), var.S##y(), var.S##z());    \
    dvar.add_##S(dx, dy, dz);                                                                                     \
    if (C == deS::bc) {                                                                                           \
      dvar.sub_L1(dx, dy, dz);                                                                                    \
      auto [nex, ney, nez] = deSitter_e_vec(var.S##x(), var.S##y(), var.S##z(), var.L1x(), var.L1y(), var.L1z()); \
      dvar.add_e1(cross_with_coef(d.S##L1_Omega(), nex, ney, nez, var.e1x(), var.e1y(), var.e1z()));              \
    }                                                                                                             \
  }

#define EVOLVE_SEL_OUT(C, S)                                                                                  \
  if (C != deS::off) {                                                                                        \
    auto [dx, dy, dz] =                                                                                       \
        cross_with_coef(d.S##L2_Omega(), d.L2x(), d.L2y(), d.L2z(), var.S##x(), var.S##y(), var.S##z());      \
    dvar.add_##S(dx, dy, dz);                                                                                 \
    if (C == deS::bc) {                                                                                       \
      if (!single_ave) {                                                                                      \
        dvar.sub_L2(dx, dy, dz);                                                                              \
        auto [nex, ney, nez] = deSitter_e_vec(var.S##x(), var.S##y(), var.S##z(), d.L2x(), d.L2y(), d.L2z()); \
        dvar.add_e2(cross_with_coef(d.S##L2_Omega(), nex, ney, nez, var.e2x(), var.e2y(), var.e2z()));        \
      } else {                                                                                                \
        dvar.add_v(SA_back_reaction(d.S##L2_Omega(), var.S##x(), var.S##y(), var.S##z(), var));               \
      }                                                                                                       \
    }                                                                                                         \
  }

template <typename Args, typename Container>
inline void deSitter_precession(bool single_ave, SLstat const &stat, Args const &args, Container const &var,
                                Container &dvar, double t) {
  using deArgs = deSitter_arg<Args, Container>;
  deArgs d{single_ave, stat, args, var};  // calculate the Omega and L2(Single average case)

  if (var.spin_num >= 1) {
    EVOLVE_SEL_IN(stat.Sin_Lin, S1);
    EVOLVE_SEL_OUT(stat.Sin_Lout, S1);
  }

  if (var.spin_num >= 2) {
    EVOLVE_SEL_IN(stat.Sin_Lin, S2);
    EVOLVE_SEL_OUT(stat.Sin_Lout, S2);
  }

  if (var.spin_num == 3) {
    EVOLVE_SEL_IN(stat.Sout_Lin, S3);
    EVOLVE_SEL_OUT(stat.Sout_Lout, S3);
  }

  if (stat.LL == deS::on) {
    dvar.add_L1(cross_with_coef(d.LL(), d.L2x(), d.L2y(), d.L2z(), var.L1x(), var.L1y(), var.L1z()));  // evolve L1
    dvar.add_e1(cross_with_coef(d.LL(), d.L2x(), d.L2y(), d.L2z(), var.e1x(), var.e1y(), var.e1z()));  // evolve e1
  }
}

}  // namespace secular
#endif
