#ifndef SECULAR_H
#define SECULAR_H

#include <iostream>
#include <array>
#include <tuple>
#include <memory>


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
        return std::to_string(task_id) + str_ave[DA] + str_pole[ctrl.Oct] + str_gr[ctrl.GR] + str_gw[ctrl.GW] + str_sl[ctrl.SL] + str_ll[ctrl.LL] +str_s[spin_num] + "\n";
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

/*template< typename ...Func>
inline auto serilize(Func...func) {
    return [=](auto const& args, auto const&x, auto &dxdt, double t) {
        (func(args, x, dxdt,t),...);
    };
}*/

    template<typename Args, typename ...Func>
    inline auto serilize(Args const &args, Func...func) {
        return [=, &args](auto const &x, auto &dxdt, double t) {
            (func(args, x, dxdt, t), ...);
        };
    }

    template<typename Ctrl>
    int ctrl_to_int(Ctrl const &c) {
        return static_cast<int>(c.DA)
               + (static_cast<int>(c.Oct) << 1)
               + (static_cast<int>(c.GR) << 2)
               + (static_cast<int>(c.GW) << 3)
               + (static_cast<int>(c.LL) << 4)
               + (static_cast<int>(c.SL) << 5);
    }

    template<typename Ctrl, typename Args, typename Container>
    auto Static_dispatch(Ctrl const &c, Args const &args) {
        int patch_id = ctrl_to_int(c);
        std::cout << patch_id << "!!\n";
        switch (patch_id) {
            case 0:
                return secular::serilize(args, single_aved_LK<false, Args, Container>);
                break;

            case 1:
                return secular::serilize(args, double_aved_LK<false, Args, Container>);
                break;

            default:
                std::cout << "wrong dispatch num~\n";
                exit(0);
        }
    }

} // namespace secular
#endif
