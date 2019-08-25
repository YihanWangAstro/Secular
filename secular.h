#ifndef SECULAR_H
#define SECULAR_H

#include <iostream>
#include <array>
#include <tuple>
#include <memory>
#include <sstream>

#include "LK.h"
#include "relativistic.h"
#include "deSitter.h"
#include "SpaceHub/src/dev-tools.hpp"

namespace secular {

    struct OrbitArgs {
        template<typename Iter>
        OrbitArgs(Iter iter, bool DA, size_t spin_num) {
            m1 = *iter, iter++;
            m2 = *iter, iter++;
            m3 = *iter, iter++;
            a_in = *iter, iter++;
            a_out = *iter, iter++;
            e_in = *iter, iter++;
            e_out = *iter, iter++;
            omega_in = *iter, iter++;
            omega_out = *iter, iter++;
            Omega_in = *iter, iter++;
            Omega_out = Omega_in - 180.0;
            i_in = *iter, iter++;
            i_out = *iter, iter++;
            if(!DA) {
                M_nu = *iter, iter++;
            }

            deg_to_rad(*this);

            s.reserve(spin_num);
            for(size_t i = 0 ; i < spin_num; ++i){
                s.emplace_back(Vec3d(*(iter+i*3), *(iter+i*3+1), *(iter+i*3+2)));
            }
        }

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
        double M_nu{0};
        std::vector <Vec3d> s;
    };

    struct Controler {
        template<typename Iter>
        Controler(Iter iter, bool _DA) {
            DA = _DA;
            Oct = static_cast<bool>(*iter), iter++;
            GR = static_cast<bool>(*iter), iter++;
            GW = static_cast<bool>(*iter), iter++;
            SL = static_cast<bool>(*iter), iter++;
            LL = static_cast<bool>(*iter);
        }

        bool DA;
        bool Oct;
        bool GR;
        bool GW;
        bool SL;
        bool LL;
    };

    const std::string str_ave[2]={":SA",":DA"};
    const std::string str_pole[2]={"|quad", "| oct"};
    const std::string str_gr[2]={"", "|GR"};
    const std::string str_gw[2]={"", "|GW"};
    const std::string str_sl[2]={"", "|S_{in}L_{out}"};
    const std::string str_ll[2]={"", "|LL"};
    const std::string str_s[4] ={"", "|S_{1}L_{in}", "|S_{1}L_{in}|S_{2}L_{in}", "|S_{1}L_{in}|S_{2}L_{in}|S_{3}L_{out}"};

    std::string get_log_title(size_t task_id, bool DA, Controler const& ctrl, size_t spin_num){
        return std::to_string(task_id) + str_ave[DA] + str_pole[ctrl.Oct] + str_gr[ctrl.GR] + str_gw[ctrl.GW] + str_sl[ctrl.SL] + str_ll[ctrl.LL] +str_s[spin_num];
    }

    enum class StopFlag {
        shrink, eof, input_err
    };

    inline double deSitter_coef(double m_self, double m_other) {
        return 0.5*consts::G/(consts::C * consts::C) * (4 + 3*m_other/m_self);
    };

    struct SecularConst {

        SecularConst(Controler const &ctrl, double _m1, double _m2, double _m3) : m1{_m1}, m2{_m2}, m3{_m3} {
            double const m12 = m1 + m2;

            mu[0] = m1 * m2 / m12;
            mu[1] = m12 * m3 / (m12 + m3);
            a_coef[0] = 1 / (consts::G * m12) / mu[0] / mu[0];
            a_coef[1] = 1 / (consts::G * (m12 + m3)) / mu[1] / mu[1];

            if (!ctrl.DA) {
                SA_acc_coef = consts::G * m3 / mu[1];
            }

            if (ctrl.GR) {
                GR_coef = 3 * pow(consts::G*m12, 1.5) / (consts::C * consts::C);
            }

            if (ctrl.GW) {
                constexpr double C5 =  consts::C *consts::C* consts::C* consts::C*consts::C;
                constexpr double G3 =  consts::G* consts::G*consts::G;

                GW_L_coef = -6.4 * pow(consts::G, 3.5) * mu[0] * mu[0] * pow(m12, 2.5) / C5;
                GW_e_coef = -304.0 / 15 * G3 * mu[0] * m12 * m12 / C5 ;
            }

            SL[0][0] = deSitter_coef(m1, m2);
            SL[1][0] = deSitter_coef(m2, m1);

            if(ctrl.SL) {
                SL[0][1] = deSitter_coef(m12, m3);
                SL[1][1] = SL[0][1];
            }

            SL[2][0] = 0;
            SL[2][1] = deSitter_coef(m3, m12);

            if(ctrl.LL) {
                LL = deSitter_coef(m12, m3);
            }
        }

    public:
        double m1;
        double m2;
        double m3;
        double mu[2];
        double a_coef[2];
        double SL[3][2];
        double LL{0};
        double SA_acc_coef{0};
        double GR_coef{0};
        double GW_L_coef{0};
        double GW_e_coef{0};
    };

    template<typename Container>
    struct Dynamic_dispatch {
        Dynamic_dispatch(Controler const &_ctrl, SecularConst const &_args) : ctrl{&_ctrl}, args{&_args} {}

        void operator()(Container const &x, Container &dxdt, double t) {
            if (ctrl->DA == true) {
                if (ctrl->Oct == true) {
                    double_aved_LK<true, SecularConst, Container>(*args, x, dxdt, t);
                } else {
                    double_aved_LK<false, SecularConst, Container>(*args, x, dxdt, t);
                }

                if (ctrl->LL == true) {
                    if (ctrl->SL == true) {
                        deSitter_precession<true, true, true, SecularConst, Container>(*args, x, dxdt, t);
                    } else {
                        deSitter_precession<true, true, false, SecularConst, Container>(*args, x, dxdt, t);
                    }
                } else {
                    if (ctrl->SL == true) {
                        deSitter_precession<true, false, true, SecularConst, Container>(*args, x, dxdt, t);
                    } else {
                        deSitter_precession<true, false, false, SecularConst, Container>(*args, x, dxdt, t);
                    }
                }
            } else {
                if (ctrl->Oct == true) {
                    single_aved_LK<true, SecularConst, Container>(*args, x, dxdt, t);
                } else {
                    single_aved_LK<false, SecularConst, Container>(*args, x, dxdt, t);
                }

                if (ctrl->LL == true) {
                    if (ctrl->SL == true) {
                        deSitter_precession<false, true, true, SecularConst, Container>(*args, x, dxdt, t);
                    } else {
                        deSitter_precession<false, true, false, SecularConst, Container>(*args, x, dxdt, t);
                    }
                } else {
                    if (ctrl->SL == true) {
                        deSitter_precession<false, false, true, SecularConst, Container>(*args, x, dxdt, t);
                    } else {
                        deSitter_precession<false, false, false, SecularConst, Container>(*args, x, dxdt, t);
                    }
                }
            }

            if (ctrl->GR == true) {
                GR_precession(*args, x, dxdt, t);
            }

            if (ctrl->GW == true) {
                GW_radiation(*args, x, dxdt, t);
            }
        }

        Controler const *ctrl;
        SecularConst const *args;
    };

/*
    template<typename Container, bool DA, bool Oct, bool GR, bool GW, bool SL, bool LL>
    struct Static_functor{
        Static_functor(SecularConst const &_args) : args{&_args} {}

        void operator()(Container const &x, Container &dxdt, double t) {
            if constexpr(DA == true) {
                if constexpr(Oct == true) {
                    double_aved_LK<true, SecularConst, Container>(*args, x, dxdt, t);
                } else {
                    double_aved_LK<false, SecularConst, Container>(*args, x, dxdt, t);
                }

                if constexpr (LL == true) {
                    if constexpr(SL == true) {
                        deSitter_precession<true, true, true, SecularConst, Container>(*args, x, dxdt, t);
                    } else {
                        deSitter_precession<true, true, false, SecularConst, Container>(*args, x, dxdt, t);
                    }
                } else {
                    if constexpr(SL == true) {
                        deSitter_precession<true, false, true, SecularConst, Container>(*args, x, dxdt, t);
                    } else {
                        deSitter_precession<true, false, false, SecularConst, Container>(*args, x, dxdt, t);
                    }
                }
            } else {
                if constexpr(Oct == true) {
                    single_aved_LK<true, SecularConst, Container>(*args, x, dxdt, t);
                } else {
                    single_aved_LK<false, SecularConst, Container>(*args, x, dxdt, t);
                }

                if constexpr(LL == true) {
                    if constexpr(SL == true) {
                        deSitter_precession<false, true, true, SecularConst, Container>(*args, x, dxdt, t);
                    } else {
                        deSitter_precession<false, true, false, SecularConst, Container>(*args, x, dxdt, t);
                    }
                } else {
                    if constexpr(SL == true) {
                        deSitter_precession<false, false, true, SecularConst, Container>(*args, x, dxdt, t);
                    } else {
                        deSitter_precession<false, false, false, SecularConst, Container>(*args, x, dxdt, t);
                    }
                }
            }

            if constexpr(GR == true) {
                GR_precession(*args, x, dxdt, t);
            }

            if constexpr(GW == true) {
                GW_radiation(*args, x, dxdt, t);
            }
        }

        SecularConst const *args;
    };

    int ctrl_to_int(Controler  const &c) {
        return static_cast<int>(c.DA)
               + (static_cast<int>(c.Oct) << 1)
               + (static_cast<int>(c.GR) << 2)
               + (static_cast<int>(c.GW) << 3)
               + (static_cast<int>(c.SL) << 4)
               + (static_cast<int>(c.LL) << 5);
    }


#define STATIC_DISPATH(ctrl, args, expr)       {                                                                    \
        int f_id = ctrl_to_int(ctrl);                                                                             \
                                                                                                                  \
        switch (f_id) {                                                                                           \
          case 0:                                                                                                 \
              {auto func = std::move(secular::Static_functor<Container, false, false, false, false, false, false>(args)); expr; break;}\
          case 1:                                                                                                 \
              {auto func = std::move(secular::Static_functor<Container, true, false, false, false, false, false>(args)); expr;break;}\
          case 2:                                                                                                 \
              {auto func = std::move(secular::Static_functor<Container, false, true, false, false, false, false>(args)); expr;break;}\
          case 3:                                                                                                 \
              {auto func = std::move(secular::Static_functor<Container, true, true, false, false, false, false>(args)); expr;break;}\
          case 4:                                                                                                 \
              {auto func = std::move(secular::Static_functor<Container, false, false, true, false, false, false>(args)); expr;break;}\
          case 5:                                                                                                 \
              {auto func = std::move(secular::Static_functor<Container, true, false, true, false, false, false>(args)); expr;break;}\
          case 6:                                                                                                 \
              {auto func = std::move(secular::Static_functor<Container, false, true, true, false, false, false>(args)); expr;break;}\
          case 7:                                                                                                 \
              {auto func = std::move(secular::Static_functor<Container, true, true, true, false, false, false>(args)); expr;break;}\
          case 8:                                                                                                 \
              {auto func = std::move(secular::Static_functor<Container, false, false, false, true, false, false>(args)); expr;break;}\
          case 9:                                                                                                 \
              {auto func = std::move(secular::Static_functor<Container, true, false, false, true, false, false>(args)); expr;break;}\
          case 10:                                                                                                 \
              {auto func = std::move(secular::Static_functor<Container, false, true, false, true, false, false>(args)); expr;break;}\
          case 11:                                                                                                 \
              {auto func = std::move(secular::Static_functor<Container, true, true, false, true, false, false>(args)); expr;break;}\
          case 12:                                                                                                 \
              {auto func = std::move(secular::Static_functor<Container, false, false, true, true, false, false>(args)); expr;break;}\
          case 13:                                                                                                 \
              {auto func = std::move(secular::Static_functor<Container, true, false, true, true, false, false>(args)); expr;break;}\
          case 14:                                                                                                 \
              {auto func = std::move(secular::Static_functor<Container, false, true, true, true, false, false>(args)); expr;break;}\
          case 15:                                                                                                 \
              {auto func = std::move(secular::Static_functor<Container, true, true, true, true, false, false>(args)); expr;break;}\
          case 16:                                                                                                 \
              {auto func = std::move(secular::Static_functor<Container, false, false, false, false, true, false>(args)); expr;break;}\
          case 17:                                                                                                 \
              {auto func = std::move(secular::Static_functor<Container, true, false, false, false, true, false>(args)); expr;break;}\
          case 18:                                                                                                 \
              {auto func = std::move(secular::Static_functor<Container, false, true, false, false, true, false>(args)); expr;break;}\
          case 19:                                                                                                 \
              {auto func = std::move(secular::Static_functor<Container, true, true, false, false, true, false>(args)); expr;break;}\
          case 20:                                                                                                 \
              {auto func = std::move(secular::Static_functor<Container, false, false, true, false, true, false>(args)); expr;break;}\
          case 21:                                                                                                 \
              {auto func = std::move(secular::Static_functor<Container, true, false, true, false, true, false>(args)); expr;break;}\
          case 22:                                                                                                 \
              {auto func = std::move(secular::Static_functor<Container, false, true, true, false, true, false>(args)); expr;break;}\
          case 23:                                                                                                 \
              {auto func = std::move(secular::Static_functor<Container, true, true, true, false, true, false>(args)); expr;break;}\
          case 24:                                                                                                 \
              {auto func = std::move(secular::Static_functor<Container, false, false, false, true, true, false>(args)); expr;break;}\
          case 25:                                                                                                 \
              {auto func = std::move(secular::Static_functor<Container, true, false, false, true, true, false>(args)); expr;break;}\
          case 26:                                                                                                 \
              {auto func = std::move(secular::Static_functor<Container, false, true, false, true, true, false>(args)); expr;break;}\
          case 27:                                                                                                 \
              {auto func = std::move(secular::Static_functor<Container, true, true, false, true, true, false>(args)); expr;break;}\
          case 28:                                                                                                 \
              {auto func = std::move(secular::Static_functor<Container, false, false, true, true, true, false>(args)); expr;break;}\
          case 29:                                                                                                 \
              {auto func = std::move(secular::Static_functor<Container, true, false, true, true, true, false>(args)); expr;break;}\
          case 30:                                                                                                 \
              {auto func = std::move(secular::Static_functor<Container, false, true, true, true, true, false>(args)); expr;break;}\
          case 31:                                                                                                 \
              {auto func = std::move(secular::Static_functor<Container, true, true, true, true, true, false>(args)); expr;break;}\
              case 32:                                                                                                 \
                  {auto func = std::move(secular::Static_functor<Container, false, false, false, false, false, true>(args)); expr;break;}\
              case 33:                                                                                                 \
                  {auto func = std::move(secular::Static_functor<Container, true, false, false, false, false, true>(args)); expr;break;}\
              case 34:                                                                                                 \
                  {auto func = std::move(secular::Static_functor<Container, false, true, false, false, false, true>(args)); expr;break;}\
              case 35:                                                                                                 \
                  {auto func = std::move(secular::Static_functor<Container, true, true, false, false, false, true>(args)); expr;break;}\
              case 36:                                                                                                 \
                  {auto func = std::move(secular::Static_functor<Container, false, false, true, false, false, true>(args)); expr;break;}\
              case 37:                                                                                                 \
                  {auto func = std::move(secular::Static_functor<Container, true, false, true, false, false, true>(args)); expr;break;}\
              case 38:                                                                                                 \
                  {auto func = std::move(secular::Static_functor<Container, false, true, true, false, false, true>(args)); expr;break;}\
              case 39:                                                                                                 \
                  {auto func = std::move(secular::Static_functor<Container, true, true, true, false, false, true>(args)); expr;break;}\
              case 40:                                                                                                 \
                  {auto func = std::move(secular::Static_functor<Container, false, false, false, true, false, true>(args)); expr;break;}\
              case 41:                                                                                                 \
                  {auto func = std::move(secular::Static_functor<Container, true, false, false, true, false, true>(args)); expr;break;}\
              case 42:                                                                                                 \
                  {auto func = std::move(secular::Static_functor<Container, false, true, false, true, false, true>(args)); expr;break;}\
              case 43:                                                                                                 \
                  {auto func = std::move(secular::Static_functor<Container, true, true, false, true, false, true>(args)); expr;break;}\
              case 44:                                                                                                 \
                  {auto func = std::move(secular::Static_functor<Container, false, false, true, true, false, true>(args)); expr;break;}\
              case 45:                                                                                                 \
                  {auto func = std::move(secular::Static_functor<Container, true, false, true, true, false, true>(args));expr; break;}\
              case 46:                                                                                                 \
                  {auto func = std::move(secular::Static_functor<Container, false, true, true, true, false, true>(args)); expr;break;}\
              case 47:                                                                                                 \
                  {auto func = std::move(secular::Static_functor<Container, true, true, true, true, false, true>(args)); expr;break;}\
              case 48:                                                                                                 \
                  {auto func = std::move(secular::Static_functor<Container, false, false, false, false, true, true>(args)); expr;break;}\
              case 49:                                                                                                 \
                  {auto func = std::move(secular::Static_functor<Container, true, false, false, false, true, true>(args)); expr;break;}\
              case 50:                                                                                                 \
                  {auto func = std::move(secular::Static_functor<Container, false, true, false, false, true, true>(args)); expr;break;}\
              case 51:                                                                                                 \
                  {auto func = std::move(secular::Static_functor<Container, true, true, false, false, true, true>(args)); expr;break;}\
              case 52:                                                                                                 \
                  {auto func = std::move(secular::Static_functor<Container, false, false, true, false, true, true>(args)); expr;break;}\
              case 53:                                                                                                 \
                  {auto func = std::move(secular::Static_functor<Container, true, false, true, false, true, true>(args)); expr;break;}\
              case 54:                                                                                                 \
                  {auto func = std::move(secular::Static_functor<Container, false, true, true, false, true, true>(args)); expr;break;}\
              case 55:                                                                                                 \
                  {auto func = std::move(secular::Static_functor<Container, true, true, true, false, true, true>(args)); expr;break;}\
              case 56:                                                                                                 \
                  {auto func = std::move(secular::Static_functor<Container, false, false, false, true, true, true>(args)); expr;break;}\
              case 57:                                                                                                 \
                  {auto func = std::move(secular::Static_functor<Container, true, false, false, true, true, true>(args)); expr;break;}\
              case 58:                                                                                                 \
                  {auto func = std::move(secular::Static_functor<Container, false, true, false, true, true, true>(args)); expr;break;}\
              case 59:                                                                                                 \
                  {auto func = std::move(secular::Static_functor<Container, true, true, false, true, true, true>(args)); expr;break;}\
              case 60:                                                                                                 \
                  {auto func = std::move(secular::Static_functor<Container, false, false, true, true, true, true>(args)); expr;break;}\
              case 61:                                                                                                 \
                  {auto func = std::move(secular::Static_functor<Container, true, false, true, true, true, true>(args)); expr;break;}\
              case 62:                                                                                                 \
                  {auto func = std::move(secular::Static_functor<Container, false, true, true, true, true, true>(args)); expr;break;}\
              case 63:                                                                                                 \
                  {auto func = std::move(secular::Static_functor<Container, true, true, true, true, true, true>(args)); expr;break;}\
        }           }                                                                                         \
*/
} // namespace secular
#endif
