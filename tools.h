#ifndef SECULAR_TOOL_H
#define SECULAR_TOOL_H

#include <algorithm>
#include <cmath>
#include <sstream>
#include <tuple>
#include <vector>
namespace secular {

#define STD_ACCESSOR(TYPE, NAME, MEMBER)  \
  inline TYPE &NAME() { return MEMBER; }; \
  inline TYPE const &NAME() const { return MEMBER; };

#define READ_GETTER(TYPE, NAME, MEMBER) \
  inline TYPE const &NAME() const { return MEMBER; };

#define OPT_STD_GETTER(COND, TYPE, NAME, MEMBER)   \
  inline TYPE &NAME() {                            \
    static_assert(COND, "method is not defined!"); \
    return MEMBER;                                 \
  };                                               \
  inline TYPE const &NAME() const {                \
    static_assert(COND, "method is not defined!"); \
    return MEMBER;                                 \
  };

#define OPT_READ_GETTER(COND, TYPE, NAME, MEMBER)  \
  inline TYPE const &NAME() const {                \
    static_assert(COND, "method is not defined!"); \
    return MEMBER;                                 \
  };

#define STD_3WAY_SETTER(NAME, X, Y, Z)                                             \
  inline void set_##NAME(double x, double y, double z) { X = x, Y = y, Z = z; }    \
  inline void add_##NAME(double x, double y, double z) { X += x, Y += y, Z += z; } \
  inline void sub_##NAME(double x, double y, double z) { X -= x, Y -= y, Z -= z; } \
  inline void set_##NAME(std::tuple<double, double, double> const &t) {            \
    X = std::get<0>(t), Y = std::get<1>(t), Z = std::get<2>(t);                    \
  }                                                                                \
  inline void add_##NAME(std::tuple<double, double, double> const &t) {            \
    X += std::get<0>(t), Y += std::get<1>(t), Z += std::get<2>(t);                 \
  }                                                                                \
  inline void sub_##NAME(std::tuple<double, double, double> const &t) {            \
    X -= std::get<0>(t), Y -= std::get<1>(t), Z -= std::get<2>(t);                 \
  }

#define OPT_3WAY_SETTER(COND, NAME, X, Y, Z)                            \
  inline void set_##NAME(double x, double y, double z) {                \
    static_assert(COND, "method is not defined!");                      \
    X = x, Y = y, Z = z;                                                \
  }                                                                     \
  inline void add_##NAME(double x, double y, double z) {                \
    static_assert(COND, "method is not defined!");                      \
    X += x, Y += y, Z += z;                                             \
  }                                                                     \
  inline void sub_##NAME(double x, double y, double z) {                \
    static_assert(COND, "method is not defined!");                      \
    X -= x, Y -= y, Z -= z;                                             \
  }                                                                     \
  inline void set_##NAME(std::tuple<double, double, double> const &t) { \
    static_assert(COND, "method is not defined!");                      \
    X = std::get<0>(t), Y = std::get<1>(t), Z = std::get<2>(t);         \
  }                                                                     \
  inline void add_##NAME(std::tuple<double, double, double> const &t) { \
    static_assert(COND, "method is not defined!");                      \
    X += std::get<0>(t), Y += std::get<1>(t), Z += std::get<2>(t);      \
  }                                                                     \
  inline void sub_##NAME(std::tuple<double, double, double> const &t) { \
    static_assert(COND, "method is not defined!");                      \
    X -= std::get<0>(t), Y -= std::get<1>(t), Z -= std::get<2>(t);      \
  }

template <typename Iter, size_t... I>
inline auto _unpack_array_(Iter iter, std::index_sequence<I...>) {
  return std::make_tuple(*(iter + I)...);
}

template <size_t num, typename Iter>
inline auto unpack_args(Iter iter) {
  return _unpack_array_(iter, std::make_index_sequence<num>());
}

template <typename Iter, typename Tup, size_t... I>
inline auto _cast_unpack_(Iter iter, Tup const &t, std::index_sequence<I...>) {
  return std::make_tuple(static_cast<decltype(std::get<I>(t))>(*(iter + I))...);
}

template <typename Iter, typename... Args>
inline auto cast_unpack(Iter iter) {
  return _cast_unpack_(iter, std::make_tuple(static_cast<Args>(0)...), std::make_index_sequence<sizeof...(Args)>());
  // return std::make_tuple(static_cast<Args>(*iter++)...);
}

void unpack_args_from_str(std::string const &str, std::vector<double> &vec, size_t num) {
  std::stringstream is{str};
  double tmp;

  vec.reserve(num);
  for (size_t i = 0; i < num; ++i) {
    is >> tmp;
    vec.emplace_back(tmp);
  }
}

#define UNPACK3(x) std::get<0>(x), std::get<1>(x), std::get<2>(x)

namespace consts {
constexpr double pi = 3.14159265358979323;
constexpr double G = 4 * pi * pi;
constexpr double C = 6.32397263e4;
constexpr double r_G_sqrt = 1.0 / (2 * pi);

constexpr double year = 1;
}  // namespace consts

enum class ReturnFlag { input_err, max_iter, finish };

bool case_insens_equals(std::string const &a, std::string const &b) {
  return std::equal(a.begin(), a.end(), b.begin(), b.end(), [](char a, char b) { return tolower(a) == tolower(b); });
}

bool is_number(const std::string &s) {
  return !s.empty() && std::find_if(s.begin(), s.end(), [](char c) { return !std::isdigit(c); }) == s.end();
}

bool str_to_bool(std::string const &key) {
  if (case_insens_equals(key, "true") || case_insens_equals(key, "1") || case_insens_equals(key, "on")) {
    return true;
  } else if (case_insens_equals(key, "false") || case_insens_equals(key, "0") || case_insens_equals(key, "off")) {
    return false;
  } else {
    throw ReturnFlag::input_err;
  }
}

bool is_on(double x) { return x > 5e-15; }

template <typename Container>
struct spin_num {
  static constexpr size_t size{Container::s_num};
};

using Tup3d = std::tuple<double, double, double>;

inline double norm2(double x, double y, double z) { return x * x + y * y + z * z; }

inline double norm2(Tup3d const &tup) { return norm2(UNPACK3(tup)); }

inline double norm(double x, double y, double z) { return sqrt(norm2(x, y, z)); }

inline double norm(Tup3d const &tup) { return sqrt(norm2(tup)); }

inline double angle(double x, double y, double z, double i, double j, double k) {
  return acos((x * i + y * j + z * k) / (norm(x, y, z) * norm(i, j, k)));
}

inline double dot(double x1, double y1, double z1, double x2, double y2, double z2) {
  return x1 * x2 + y1 * y2 + z1 * z2;
}

inline double dot(Tup3d const &t1, Tup3d const &t2) { return dot(UNPACK3(t1), UNPACK3(t2)); }

inline auto cross(double x1, double y1, double z1, double x2, double y2, double z2) {
  return std::make_tuple(y1 * z2 - y2 * z1, z1 * x2 - z2 * x1, x1 * y2 - x2 * y1);
}

inline auto cross(Tup3d const &t1, Tup3d const &t2) { return cross(UNPACK3(t1), UNPACK3(t2)); }

inline auto cross_with_coef(double A, double x1, double y1, double z1, double x2, double y2, double z2) {
  return std::make_tuple(A * (y1 * z2 - y2 * z1), A * (z1 * x2 - z2 * x1), A * (x1 * y2 - x2 * y1));
}

inline auto cross_with_coef(double A, Tup3d const &t1, Tup3d const &t2) {
  return cross_with_coef(A, UNPACK3(t1), UNPACK3(t2));
}

double calc_angular_mom(double m_in, double m_out, double a) {
  double mu = m_in * m_out / (m_in + m_out);
  return mu * sqrt(consts::G * (m_in + m_out) * a);
}

auto unit_j(double i, double Omega) {
  double sini = sin(i);
  double cosi = cos(i);
  double sinO = sin(Omega);
  double cosO = cos(Omega);
  return std::make_tuple(sini * sinO, -sini * cosO, cosi);
}

auto unit_e(double i, double omega, double Omega) {
  double sini = sin(i);
  double cosi = cos(i);
  double sinO = sin(Omega);
  double cosO = cos(Omega);
  double sino = sin(omega);
  double coso = cos(omega);

  return std::make_tuple(coso * cosO - sino * cosi * sinO, coso * sinO + sino * cosi * cosO, sino * sini);
}

auto unit_peri_v(double i, double omega, double Omega) {
  auto [ex, ey, ez] = unit_e(i, omega, Omega);
  auto [jx, jy, jz] = unit_j(i, Omega);
  return cross(jx, jy, jz, ex, ey, ez);
}

template <typename... Args>
void deg_to_rad(Args &... args) {
  constexpr double rad = consts::pi / 180.0;
  ((args *= rad), ...);
}

inline auto calc_orbit_args(double Coef, double lx, double ly, double lz, double ex, double ey, double ez) {
  double e_sqr = norm2(ex, ey, ez);

  double j_sqr = fabs(1 - e_sqr);

  double j = sqrt(j_sqr);

  double L_norm = norm(lx, ly, lz);

  double L = L_norm / j;

  double a = Coef * L * L;

  return std::make_tuple(e_sqr, j_sqr, j, L_norm, L, a);
}

inline auto calc_a_eff(double Coef, double lx, double ly, double lz, double ex, double ey, double ez) {
  double e_sqr = norm2(ex, ey, ez);

  double j_sqr = fabs(1 - e_sqr);

  double j = sqrt(j_sqr);

  double L_sqr = norm2(lx, ly, lz);

  return Coef * L_sqr / j;
}

inline auto calc_a_j(double Coef, double lx, double ly, double lz, double ex, double ey, double ez) {
  double e_sqr = norm2(ex, ey, ez);

  double j_sqr = fabs(1 - e_sqr);

  double j = sqrt(j_sqr);

  double L_sqr = norm2(lx, ly, lz);

  return std::make_tuple(Coef * L_sqr / j_sqr, j);
}

inline auto calc_a(double Coef, double lx, double ly, double lz, double ex, double ey, double ez) {
  double e_sqr = norm2(ex, ey, ez);

  double j_sqr = fabs(1 - e_sqr);

  double L_sqr = norm2(lx, ly, lz);

  return Coef * L_sqr / j_sqr;
}

const std::string str_ave[2] = {":SA", ":DA"};
const std::string str_pole[2] = {"|quad", "| oct"};
const std::string str_gr_in[2] = {"", "|GR_{in}"};
const std::string str_gr_out[2] = {"", "|GR_{out}"};
const std::string str_gw_in[2] = {"", "|GW_{in}"};
const std::string str_gw_out[2] = {"", "|GW_{out}"};

const std::string str_sin_lin[4] = {"", "|^{dS}S_{in}L_{in}", "|^{back}S_{in}L_{in}", "|^{dS+back}S_{in}L_{in}"};
const std::string str_sin_lout[4] = {"", "|^{dS}S_{in}L_{out}", "|^{back}S_{in}L_{out}", "|^{dS+back}S_{in}L_{out}"};
const std::string str_sout_lin[4] = {"", "|^{LT}S_{out}L_{in}", "|^{back}S_{out}L_{in}", "|^{LT+back}S_{out}L_{in}"};
const std::string str_sout_lout[4] = {"", "|^{ds}S_{out}L_{out}", "|^{back}S_{out}L_{out}",
                                      "|^{ds+back}S_{out}L_{out}"};

const std::string str_sin_sin[4] = {"", "|^{LT}S_{in}S_{in}", "|^{back}S_{in}S_{in}", "|^{LT+back}S_{in}S_{in}"};
const std::string str_sin_sout[4] = {"", "|^{LT}S_{in}S_{out}", "|^{back}S_{in}S_{out}", "|^{LT+back}S_{in}S_{out}"};

const std::string str_ll[4] = {
    "",
    "|^{dS}L_{in}L_{out}",
    "|^{back}L_{in}L_{out}",
    "|^{dS+back}L_{in}L_{out}",
};

template <typename Controler>
std::string get_log_title(Controler const &ctrl) {
  return std::string{"config:"} + str_ave[to_index(ctrl.ave_method)] + str_pole[ctrl.Oct] + str_gr_in[ctrl.GR_in] +
         str_gr_out[ctrl.GR_out] + str_gw_in[ctrl.GW_in] + str_gw_out[ctrl.GW_out] +
         str_sin_lin[to_index(ctrl.Sin_Lin)] + str_sin_lout[to_index(ctrl.Sin_Lout)] +
         str_sout_lin[to_index(ctrl.Sout_Lin)] + str_sout_lout[to_index(ctrl.Sout_Lout)] +
         str_sin_sin[to_index(ctrl.Sin_Sin)] + str_sin_sout[to_index(ctrl.Sin_Sout)] + str_ll[to_index(ctrl.LL)];
}
}  // namespace secular
#endif
