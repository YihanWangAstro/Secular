#ifndef SECULAR_H
#define SECULAR_H

#include <array>
#include <memory>
#include <tuple>

#include "LK.h"
#include "SpaceHub/src/tools/config-reader.hpp"
#include "deSitter.h"
#include "relativistic.h"
#include "stellar.h"
#include "tools.h"

namespace secular {

class SecularArray : public std::array<double, 21> {
 public:
  SecularArray() = default;

  READ_GETTER(double, L1x, (*this)[0]);

  READ_GETTER(double, L1y, (*this)[1]);

  READ_GETTER(double, L1z, (*this)[2]);

  READ_GETTER(double, e1x, (*this)[3]);

  READ_GETTER(double, e1y, (*this)[4]);

  READ_GETTER(double, e1z, (*this)[5]);

  READ_GETTER(double, L2x, (*this)[6]);

  READ_GETTER(double, L2y, (*this)[7]);

  READ_GETTER(double, L2z, (*this)[8]);

  READ_GETTER(double, e2x, (*this)[9]);

  READ_GETTER(double, e2y, (*this)[10]);

  READ_GETTER(double, e2z, (*this)[11]);

  READ_GETTER(double, rx, (*this)[6]);

  READ_GETTER(double, ry, (*this)[7]);

  READ_GETTER(double, rz, (*this)[8]);

  READ_GETTER(double, vx, (*this)[9]);

  READ_GETTER(double, vy, (*this)[10]);

  READ_GETTER(double, vz, (*this)[11]);

  // READ_GETTER(auto, L1, std::tie((*this)[0], (*this)[1], (*this)[2]));

  // READ_GETTER(auto, e1, std::tie((*this)[3], (*this)[4], (*this)[5]));

  // READ_GETTER(auto, L2, std::tie((*this)[6], (*this)[7], (*this)[8]));

  // READ_GETTER(auto, e2, std::tie((*this)[9], (*this)[10], (*this)[11]));

  // READ_GETTER(auto, r, std::tie((*this)[6], (*this)[7], (*this)[8]));

  // READ_GETTER(auto, v, std::tie((*this)[9], (*this)[10], (*this)[11]));

  STD_3WAY_SETTER(L1, (*this)[0], (*this)[1], (*this)[2]);

  STD_3WAY_SETTER(e1, (*this)[3], (*this)[4], (*this)[5]);

  STD_3WAY_SETTER(L2, (*this)[6], (*this)[7], (*this)[8]);

  STD_3WAY_SETTER(e2, (*this)[9], (*this)[10], (*this)[11]);

  STD_3WAY_SETTER(r, (*this)[6], (*this)[7], (*this)[8]);

  STD_3WAY_SETTER(v, (*this)[9], (*this)[10], (*this)[11]);

  READ_GETTER(double, S1x, (*this)[12]);

  READ_GETTER(double, S1y, (*this)[13]);

  READ_GETTER(double, S1z, (*this)[14]);

  READ_GETTER(double, S2x, (*this)[15]);

  READ_GETTER(double, S2y, (*this)[16]);

  READ_GETTER(double, S2z, (*this)[17]);

  READ_GETTER(double, S3x, (*this)[18]);

  READ_GETTER(double, S3y, (*this)[19]);

  READ_GETTER(double, S3z, (*this)[20]);

  STD_3WAY_SETTER(S1, (*this)[12], (*this)[13], (*this)[14]);

  STD_3WAY_SETTER(S2, (*this)[15], (*this)[16], (*this)[17]);

  STD_3WAY_SETTER(S3, (*this)[18], (*this)[19], (*this)[20]);

  friend std::ostream &operator<<(std::ostream &os, SecularArray const &arr) {
    for (auto a : arr) {
      os << a << ' ';
    }
    return os;
  }

  auto spin_begin() { return this->begin() + 12; }
};

struct Controller {
  double stop_a_in() const { return GW_stop_a_; }

  double GW_in_ratio;
  double GW_out_ratio;
  LK_method ave_method;
  bool Quad{true};
  bool Oct{false};
  bool GR_in;
  bool GR_out;
  bool GW_in;
  bool GW_out;
  deS Sin_Lin;
  deS Sin_Lout;
  deS Sout_Lin;
  deS Sout_Lout;
  deS Sin_Sin;
  deS Sin_Sout;
  deS LL;

  void set_stop_a_in(double a_stop) { GW_stop_a_ = a_stop; }

  Controller(space::tools::ConfigReader &cfg) {
    ave_method = str_to_LK_enum(cfg.get<std::string>("LK_method"));

    Quad = str_to_bool(cfg.get<std::string>("quad"));

    Oct = str_to_bool(cfg.get<std::string>("oct"));

    GR_in = str_to_bool(cfg.get<std::string>("GR_in"));

    GR_out = str_to_bool(cfg.get<std::string>("GR_out"));

    GW_in_ratio = cfg.get<double>("GW_in");

    GW_in = is_on(GW_in_ratio);

    GW_out_ratio = cfg.get<double>("GW_out");

    GW_out = is_on(GW_out_ratio);

    Sin_Lin = str_to_spin_orbit_enum(cfg.get<std::string>("Sin_Lin"));

    Sin_Lout = str_to_spin_orbit_enum(cfg.get<std::string>("Sin_Lout"));

    Sout_Lin = str_to_spin_orbit_enum(cfg.get<std::string>("Sout_Lin"));

    Sout_Lout = str_to_spin_orbit_enum(cfg.get<std::string>("Sout_Lout"));

    Sin_Sin = str_to_spin_orbit_enum(cfg.get<std::string>("Sin_Sin"));

    Sin_Sout = str_to_spin_orbit_enum(cfg.get<std::string>("Sin_Sout"));

    LL = str_to_spin_orbit_enum(cfg.get<std::string>("LL"));
  }

  std::string initial_format() {
    std::string base =
        "task_id  t_{end}[yr]  dt_{output}[yr]  m_{1}[m_{solar}]  m_{2}[m_{solar}]  m_{3}[m_{solar}]  a_{in}[au]  "
        "a_{out}[au]  e_{in}  e_{out}  omega_{in}[deg]  omega_{out}[deg]  Omega[deg]  i_{in}[deg]  i_{out}[deg]";

    if (need_anomaly()) {
      base += "  M(mean anomaly)[deg]";
    }

    if (need_S_in()) {
      base += "  S_{1,x}  S_{1,y}  S_{1,z}  S_{2,x}  S_{2,y}  S_{2,z}";
    }

    if (need_S_out()) {
      base += "  S_{3,x}  S_{3,y}  S_{3,z}";
    }

    return base;
  }

  std::string output_format() {
    std::string base = "t[yr]  L_{1,x}[?]  L_{1,y}[?]  L_{1,z}[?]  e_{1,x}  e_{1,y}  e_{1,z}";

    if (need_anomaly()) {
      base += "  r_{out,x}[?]  r_{out,y}[?]  r_{out,z}[?]  v_{out,x}[?]  v_{out,y}[?]  v_{out,z}[?]";
    } else {
      base += "  L_{2,x}[?]  L_{2,y}[?]  L_{2,z}[?]  e_{2,x}  e_{2,y}  e_{2,z}";
    }

    if (need_S_in()) {
      base += "  S_{1,x}  S_{1,y}  S_{1,z}  S_{2,x}  S_{2,y}  S_{2,z}";
    }

    if (need_S_out()) {
      base += "  S_{3,x}  S_{3,y}  S_{3,z}";
    }

    return base;
  }

  bool need_anomaly() { return ave_method == LK_method::SA; }

  bool need_S_in() {
    return (Sin_Lin != deS::off) || (Sin_Lout != deS::off) || (Sin_Sin != deS::off) || (Sin_Sout != deS::off);
  }

  bool need_S_out() { return (Sout_Lin != deS::off) || (Sout_Lout != deS::off) || (Sin_Sout != deS::off); }

 private:
  double GW_stop_a_;
};

class SecularConst {
 public:
  SecularConst(double _m1, double _m2, double _m3) : stellar_{_m1, _m2, _m3} { calculate_coef(_m1, _m2, _m3); }

  SecularConst() = default;

  READ_GETTER(double, m1, stellar_.m1());

  READ_GETTER(double, m2, stellar_.m2());

  READ_GETTER(double, m3, stellar_.m3());

  READ_GETTER(double, m1_age, stellar_.m1_age());

  READ_GETTER(double, m2_age, stellar_.m2_age());

  READ_GETTER(double, m3_age, stellar_.m3_age());

  READ_GETTER(bool, m1_dead, stellar_.m1_dead());

  READ_GETTER(bool, m2_dead, stellar_.m2_dead());

  READ_GETTER(bool, m3_dead, stellar_.m3_dead());

  READ_GETTER(double, m12, basic_.m12());

  READ_GETTER(double, m_tot, basic_.m_tot());

  READ_GETTER(double, mu_in, basic_.mu_in());

  READ_GETTER(double, mu_out, basic_.mu_out());

  READ_GETTER(double, a_in_coef, basic_.a_in_coef());

  READ_GETTER(double, a_out_coef, basic_.a_out_coef());

  READ_GETTER(double, SA_acc_coef, basic_.SA_acc_coef());

  READ_GETTER(double, GR_in_coef, GR_in_.GR_coef());

  READ_GETTER(double, GW_L_in_coef, GR_in_.GW_L_coef());

  READ_GETTER(double, GW_e_in_coef, GR_in_.GW_L_coef());

  READ_GETTER(double, GR_out_coef, GR_out_.GR_coef());

  READ_GETTER(double, GW_L_out_coef, GR_out_.GW_L_coef());

  READ_GETTER(double, GW_e_out_coef, GR_out_.GW_L_coef());

  READ_GETTER(double, LL, SL_.LL());

  READ_GETTER(double, S1L1, SL_.S1L1());

  READ_GETTER(double, S1L2, SL_.S1L2());

  READ_GETTER(double, S2L1, SL_.S2L1());

  READ_GETTER(double, S2L2, SL_.S2L2());

  READ_GETTER(double, S3L1, SL_.S3L1());

  READ_GETTER(double, S3L2, SL_.S3L2());

  READ_GETTER(double, S1S2, SL_.S1S2());

  READ_GETTER(double, S1S3, SL_.S1S3());

  READ_GETTER(double, S2S3, SL_.S2S3());

  void make_m1_exploded() {
    stellar_.make_m1_exploded();
    this->calculate_coef(m1(), m2(), m3());
  }

  void make_m2_exploded() {
    stellar_.make_m2_exploded();
    this->calculate_coef(m1(), m2(), m3());
  }

  void make_m3_exploded() {
    stellar_.make_m3_exploded();
    this->calculate_coef(m1(), m2(), m3());
  }

 private:
  BasicConst basic_;
  GRConst GR_in_;
  GRConst GR_out_;
  SLConst SL_;
  StellarConst stellar_;

  void calculate_coef(double _m1, double _m2, double _m3) {
    basic_.calculate_coef(_m1, _m2, _m3);
    GR_in_.calculate_coef(basic_.m12(), basic_.mu_in());
    GR_out_.calculate_coef(basic_.m_tot(), basic_.mu_out());
    SL_.calculate_coef(_m1, _m2, _m3);
  };
};

template <typename Container>
struct Dynamic_dispatch {
  using ConstArg = SecularConst;

  Dynamic_dispatch(Controller const &_ctrl, ConstArg const &_args) : ctrl{&_ctrl}, args{&_args} {}

  void operator()(Container const &x, Container &dxdt, double t) {
    std::fill(dxdt.begin(), dxdt.end(), 0);

    Lidov_Kozai(*ctrl, *args, x, dxdt);

    spin_orbit_coupling(*ctrl, *args, x, dxdt);

    GR_precession(*ctrl, *args, x, dxdt);

    GW_radiation(*ctrl, *args, x, dxdt);
  }

  Controller const *ctrl;
  SecularConst const *args;
};
}  // namespace secular
#endif
