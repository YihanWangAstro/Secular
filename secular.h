#ifndef SECULAR_H
#define SECULAR_H

#include <array>
#include <tuple>
#include <memory>

#include "LK.h"
#include "relativistic.h"
#include "deSitter.h"


namespace secular {

    template<size_t spin_num>
    class SecularArray : public std::array<double, 12 + 3 * spin_num> {
        static_assert(spin_num <= 3, "The max spin number is 3!");
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

        READ_GETTER(auto, L1, std::tie((*this)[0], (*this)[1], (*this)[2]));

        READ_GETTER(auto, e1, std::tie((*this)[3], (*this)[4], (*this)[5]));

        READ_GETTER(auto, L2, std::tie((*this)[6], (*this)[7], (*this)[8]));

        READ_GETTER(auto, e2, std::tie((*this)[9], (*this)[10], (*this)[11]));

        READ_GETTER(auto, r, std::tie((*this)[6], (*this)[7], (*this)[8]));

        READ_GETTER(auto, v, std::tie((*this)[9], (*this)[10], (*this)[11]));

        STD_3WAY_SETTER(L1, (*this)[0], (*this)[1], (*this)[2]);

        STD_3WAY_SETTER(e1, (*this)[3], (*this)[4], (*this)[5]);

        STD_3WAY_SETTER(L2, (*this)[6], (*this)[7], (*this)[8]);

        STD_3WAY_SETTER(e2, (*this)[9], (*this)[10], (*this)[11]);

        STD_3WAY_SETTER(r, (*this)[6], (*this)[7], (*this)[8]);

        STD_3WAY_SETTER(v, (*this)[9], (*this)[10], (*this)[11]);

        OPT_READ_GETTER(spin_num > 0, double, S1x, (*this)[12]);

        OPT_READ_GETTER(spin_num > 0, double, S1y, (*this)[13]);

        OPT_READ_GETTER(spin_num > 0, double, S1z, (*this)[14]);

        OPT_READ_GETTER(spin_num > 1, double, S2x, (*this)[15]);

        OPT_READ_GETTER(spin_num > 1, double, S2y, (*this)[16]);

        OPT_READ_GETTER(spin_num > 1, double, S2z, (*this)[17]);

        OPT_READ_GETTER(spin_num > 2, double, S3x, (*this)[18]);

        OPT_READ_GETTER(spin_num > 2, double, S3y, (*this)[19]);

        OPT_READ_GETTER(spin_num > 2, double, S3z, (*this)[20]);

        OPT_3WAY_SETTER(spin_num > 0, S1, (*this)[12], (*this)[13], (*this)[14]);

        OPT_3WAY_SETTER(spin_num > 1, S2, (*this)[15], (*this)[16], (*this)[17]);

        OPT_3WAY_SETTER(spin_num > 2, S3, (*this)[18], (*this)[19], (*this)[20]);

        auto spin_begin() {
            return this->begin() + 12;
        }

        constexpr static size_t s_num{spin_num};
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

    enum class StopFlag {
        shrink, eof, input_err
    };

    template<size_t spin_num>
    class SecularConst {
    public:
        SecularConst(double _m1, double _m2, double _m3) : basic_{_m1, _m2, _m3}, GR_{basic_.m12(), basic_.mu_in()}, SL_{_m1, _m2, _m3} {}

        SecularConst() = default;

        READ_GETTER(double, m1, basic_.m1());

        READ_GETTER(double, m2, basic_.m2());

        READ_GETTER(double, m3, basic_.m3());

        READ_GETTER(double, m12, basic_.m12());

        READ_GETTER(double, m_tot, basic_.m_tot());

        READ_GETTER(double, mu_in, basic_.mu_in());

        READ_GETTER(double, mu_out, basic_.mu_out());

        READ_GETTER(double, a_in_coef, basic_.a_in_coef());

        READ_GETTER(double, a_out_coef, basic_.a_out_coef());

        READ_GETTER(double, SA_acc_coef, basic_.SA_acc_coef());

        READ_GETTER(double, GR_coef, GR_.GR_coef());

        READ_GETTER(double, GW_L_coef, GR_.GW_L_coef());

        READ_GETTER(double, GW_e_coef, GR_.GW_L_coef());

        READ_GETTER(double, LL, SL_.LL());

        READ_GETTER(double, S1L1, SL_.S1L1());

        READ_GETTER(double, S1L2, SL_.S1L2());

        READ_GETTER(double, S2L1, SL_.S2L1());

        READ_GETTER(double, S2L2, SL_.S2L2());

        READ_GETTER(double, S3L1, SL_.S3L1());

        READ_GETTER(double, S3L2, SL_.S3L2());
    private:
        BasicConst basic_;
        GRConst GR_;
        SLConst<spin_num> SL_;
    };

    template<typename Container>
    struct Dynamic_dispatch {
        using ConstArg = SecularConst<spin_num<Container>::size>;

        Dynamic_dispatch(Controler const &_ctrl, ConstArg const &_args) : ctrl{&_ctrl}, args{&_args} {}

        void operator()(Container const &x, Container &dxdt, double t) {
            if (ctrl->DA == true) {
                if (ctrl->Oct == true) {
                    double_aved_LK<true, ConstArg, Container>(*args, x, dxdt, t);
                } else {
                    double_aved_LK<false, ConstArg, Container>(*args, x, dxdt, t);
                }

                if (ctrl->LL == true) {
                    if (ctrl->SL == true) {
                        deSitter_precession<true, SLstat<deS::on, deS::on, deS::off, deS::bc, deS::on>, ConstArg, Container>(*args, x, dxdt, t);
                    } else {
                        deSitter_precession<true, SLstat<deS::on, deS::off, deS::off, deS::bc, deS::on>, ConstArg, Container>(*args, x, dxdt, t);
                    }
                } else {
                    if (ctrl->SL == true) {
                        deSitter_precession<true, SLstat<deS::on, deS::on, deS::off, deS::bc, deS::off>, ConstArg, Container>(*args, x, dxdt, t);
                    } else {
                        deSitter_precession<true, SLstat<deS::on, deS::off, deS::off, deS::bc, deS::off>, ConstArg, Container>(*args, x, dxdt, t);
                    }
                }
            } else {
                if (ctrl->Oct == true) {
                    single_aved_LK<true, ConstArg, Container>(*args, x, dxdt, t);
                } else {
                    single_aved_LK<false, ConstArg, Container>(*args, x, dxdt, t);
                }

                if (ctrl->LL == true) {
                    if (ctrl->SL == true) {
                        deSitter_precession<false, SLstat<deS::on, deS::on, deS::off, deS::bc, deS::on>, ConstArg, Container>(*args, x, dxdt, t);
                    } else {
                        deSitter_precession<false, SLstat<deS::on, deS::off, deS::off, deS::bc, deS::on>, ConstArg, Container>(*args, x, dxdt, t);
                    }
                } else {
                    if (ctrl->SL == true) {
                        deSitter_precession<false, SLstat<deS::on, deS::on, deS::off, deS::bc, deS::off>, ConstArg, Container>(*args, x, dxdt, t);
                    } else {
                        deSitter_precession<false, SLstat<deS::on, deS::off, deS::off, deS::bc, deS::off>, ConstArg, Container>(*args, x, dxdt, t);
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
        SecularConst<spin_num<Container>::size> const *args;
    };


    template<typename Container, bool DA, bool Oct, bool GR, bool GW, bool SL, bool LL>
    struct Static_functor {
        using ConstArg = SecularConst<spin_num<Container>::size>;

        Static_functor(ConstArg const &_args) : args{&_args} {}

        void operator()(Container const &x, Container &dxdt, double t) {
            if constexpr(DA == true) {
                if constexpr(Oct == true) {
                    double_aved_LK<true, ConstArg, Container>(*args, x, dxdt, t);
                } else {
                    double_aved_LK<false, ConstArg, Container>(*args, x, dxdt, t);
                }

                if constexpr (LL == true) {
                    if constexpr(SL == true) {
                        deSitter_precession<true, SLstat<deS::on, deS::on, deS::off, deS::bc, deS::on>, ConstArg, Container>(*args, x, dxdt, t);
                    } else {
                        deSitter_precession<true, SLstat<deS::on, deS::off, deS::off, deS::bc, deS::on>, ConstArg, Container>(*args, x, dxdt, t);
                    }
                } else {
                    if constexpr(SL == true) {
                        deSitter_precession<true, SLstat<deS::on, deS::on, deS::off, deS::bc, deS::off>, ConstArg, Container>(*args, x, dxdt, t);
                    } else {
                        deSitter_precession<true, SLstat<deS::on, deS::off, deS::off, deS::bc, deS::off>, ConstArg, Container>(*args, x, dxdt, t);
                    }
                }
            } else {
                if constexpr(Oct == true) {
                    single_aved_LK<true, ConstArg, Container>(*args, x, dxdt, t);
                } else {
                    single_aved_LK<false, ConstArg, Container>(*args, x, dxdt, t);
                }

                if constexpr(LL == true) {
                    if constexpr(SL == true) {
                        deSitter_precession<false, SLstat<deS::on, deS::on, deS::off, deS::bc, deS::on>, ConstArg, Container>(*args, x, dxdt, t);
                    } else {
                        deSitter_precession<false, SLstat<deS::on, deS::off, deS::off, deS::bc, deS::on>, ConstArg, Container>(*args, x, dxdt, t);
                    }
                } else {
                    if constexpr(SL == true) {
                        deSitter_precession<false, SLstat<deS::on, deS::on, deS::off, deS::bc, deS::off>, ConstArg, Container>(*args, x, dxdt, t);
                    } else {
                        deSitter_precession<false, SLstat<deS::on, deS::off, deS::off, deS::bc, deS::off>, ConstArg, Container>(*args, x, dxdt, t);
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

        ConstArg const *args;
    };

    int ctrl_to_int(Controler const &c) {
        return static_cast<int>(c.DA)
               + (static_cast<int>(c.Oct) << 1)
               + (static_cast<int>(c.GR) << 2)
               + (static_cast<int>(c.GW) << 3)
               + (static_cast<int>(c.SL) << 4)
               + (static_cast<int>(c.LL) << 5);
    }


#define STATIC_DISPATH(ctrl, args, expr)       {                                                                  \
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

} // namespace secular
#endif
