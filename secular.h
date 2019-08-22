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

    template<size_t spin_num>
    struct OrbitArgs {
        static_assert(spin_num<=3, "Spin number is not allowed to be larger than 3!");

        static constexpr size_t SpinNum{spin_num};
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
        double M_nu;
        std::array <Vec3d, SpinNum> s;

        friend std::istream &operator>>(std::istream &is, OrbitArgs &t) {
            space::input(is, t.m1, t.m2, t.m3, t.a_in, t.a_out, t.e_in, t.e_out, t.omega_in, t.omega_out, t.Omega_in, t.i_in, t.i_out, t.M_nu);
            t.Omega_out = t.Omega_in - 180.0;
            for (auto &ss : t.s) {
                is >> ss;
            }
            return is;
        }

        friend std::ostream &operator<<(std::ostream &os, OrbitArgs const &t) {
            space::display(os, t.m1, t.m2, t.m3, t.a_in, t.a_out, t.e_in, t.e_out, t.omega_in, t.omega_out, t.Omega_in, t.i_in, t.i_out, t.M_nu);
            for (auto const &ss : t.s) {
                os << ' ' << ss;
            }
            return os;
        }
    };

    struct Controler {
        size_t id;
        bool write_traj;
        bool write_end;
        double end_time;
        double dt_out;
        bool DA;
        bool Oct;
        bool GR;
        bool GW;
        bool SL;
        bool LL;

        friend std::istream &operator>>(std::istream &is, Controler &t) {
            space::input(is, t.id, t.write_traj, t.write_end, t.end_time, t.dt_out, t.DA, t.Oct, t.GR, t.GW, t.SL, t.LL);
            return is;
        }

        friend std::ostream &operator<<(std::ostream &os, Controler const &t) {
            space::display(os, t.id, t.write_traj, t.write_end, t.end_time, t.dt_out, t.DA, t.Oct, t.GR, t.GW, t.SL, t.LL);
            return os;
        }
    };

    template<size_t spin_num>
    struct Task {
        static_assert(spin_num<=3, "Spin number is not allowed to be larger than 3!");

        Controler ctrl;
        OrbitArgs<spin_num> obt_args;

        friend std::istream &operator>>(std::istream &is, Task &t) {
            space::input(is, t.ctrl, t.obt_args);
            return is;
        }

        friend std::ostream &operator<<(std::ostream &os, Task const &t) {
            space::display(os, t.ctrl, t.obt_args);
            return os;
        }
    };

    enum class StopFlag {
        shrink
    };

    inline double deSitter_coef(double m_self, double m_other) {
        return 0.5*consts::G/(consts::C * consts::C) * (4 + 3*m_other/m_self);
    };

    template<typename Ctrl, size_t spin_num>
    struct SecularArg {
        static_assert(spin_num<=3, "Spin number is not allowed to be larger than 3!");

        SecularArg(Ctrl const &ctrl, double _m1, double _m2, double _m3) : m1{_m1}, m2{_m2}, m3{_m3} {
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

            if constexpr(spin_num >= 1)
              SL[0][0] = deSitter_coef(m1, m2);

            if constexpr(spin_num >= 2)
              SL[1][0] = deSitter_coef(m2, m1);

            if(ctrl.SL) {
                SL[0][1] = deSitter_coef(m12, m3);
                SL[1][1] = SL[0][1];
            }

            if constexpr(spin_num >= 3){
                SL[2][0] = 0;
                SL[2][1] = deSitter_coef(m3, m12);
            }

            if(ctrl.LL) {
                LL = deSitter_coef(m12, m3);
            }

            std::cout << GR_coef << "\n";
        }

    public:
        double m1;
        double m2;
        double m3;
        double mu[2];
        double a_coef[2];
        double SL[spin_num][2];
        double LL{0};
        double SA_acc_coef{0};
        double GR_coef{0};
        double GW_L_coef{0};
        double GW_e_coef{0};
    };


    template<typename Ctrl, typename Args, typename Container>
    struct Dynamic_dispatch {

        Dynamic_dispatch(Ctrl const &_ctrl, Args const &_args) : ctrl{&_ctrl}, args{&_args} {}

        void operator()(Container const &x, Container &dxdt, double t) {
            if (ctrl->DA == true) {
                if (ctrl->Oct == true) {
                    double_aved_LK<true, Args, Container>(*args, x, dxdt, t);
                } else {
                    double_aved_LK<false, Args, Container>(*args, x, dxdt, t);
                }

                if (ctrl->LL == true) {
                    if (ctrl->SL == true) {
                        deSitter_precession<true, true, true, Args, Container>(*args, x, dxdt, t);
                    } else {
                        deSitter_precession<true, true, false, Args, Container>(*args, x, dxdt, t);
                    }
                } else {
                    if (ctrl->SL == true) {
                        deSitter_precession<true, false, true, Args, Container>(*args, x, dxdt, t);
                    } else {
                        deSitter_precession<true, false, false, Args, Container>(*args, x, dxdt, t);
                    }
                }
            } else {
                if (ctrl->Oct == true) {
                    single_aved_LK<true, Args, Container>(*args, x, dxdt, t);
                } else {
                    single_aved_LK<false, Args, Container>(*args, x, dxdt, t);
                }

                if (ctrl->LL == true) {
                    if (ctrl->SL == true) {
                        deSitter_precession<false, true, true, Args, Container>(*args, x, dxdt, t);
                    } else {
                        deSitter_precession<false, true, false, Args, Container>(*args, x, dxdt, t);
                    }
                } else {
                    if (ctrl->SL == true) {
                        deSitter_precession<false, false, true, Args, Container>(*args, x, dxdt, t);
                    } else {
                        deSitter_precession<false, false, false, Args, Container>(*args, x, dxdt, t);
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

        Ctrl const *ctrl;
        Args const *args;
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
