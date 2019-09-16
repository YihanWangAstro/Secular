#ifndef SECULAR_H
#define SECULAR_H

#include <array>
#include <tuple>
#include <memory>

#include "LK.h"
#include "relativistic.h"
#include "deSitter.h"
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

        /*READ_GETTER(auto, L1, std::tie((*this)[0], (*this)[1], (*this)[2]));

        READ_GETTER(auto, e1, std::tie((*this)[3], (*this)[4], (*this)[5]));

        READ_GETTER(auto, L2, std::tie((*this)[6], (*this)[7], (*this)[8]));

        READ_GETTER(auto, e2, std::tie((*this)[9], (*this)[10], (*this)[11]));

        READ_GETTER(auto, r, std::tie((*this)[6], (*this)[7], (*this)[8]));

        READ_GETTER(auto, v, std::tie((*this)[9], (*this)[10], (*this)[11]));*/

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

        friend std::ostream& operator<<(std::ostream&os, SecularArray const& arr) {
            for(auto a : arr){
                os << a << ' ';
            }
            return os;
        }

        auto spin_begin() {
            return this->begin() + 12;
        }
    };

    struct Controler {
        double stop_a_in() const{
            return GW_stop_a_;
        }

        double GW_in_ratio;
        double GW_out_ratio;
        LK_method LK_method;
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

        void set_stop_a_in(double a_stop){
            GW_stop_a_ = a_stop;
        }
    private:
        double GW_stop_a_;
    };

    class SecularConst {
    public:
        SecularConst(double _m1, double _m2, double _m3) : basic_{_m1, _m2, _m3}, GR_in_{basic_.m12(), basic_.mu_in()}, GR_out_{basic_.m_tot(), basic_.mu_out()}, SL_{_m1, _m2, _m3} {}

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
    private:
        BasicConst basic_;
        GRConst GR_in_;
        GRConst GR_out_;
        SLConst SL_;
    };

    template<typename Container>
    struct Dynamic_dispatch {
        using ConstArg = SecularConst;

        Dynamic_dispatch(Controler const &_ctrl, ConstArg const &_args) : ctrl{&_ctrl}, args{&_args} {}

        void operator()(Container const &x, Container &dxdt, double t) {
            std::fill(dxdt.begin(), dxdt.end(), 0);

            Lidov_Kozai(*ctrl, *args, x, dxdt);

            spin_orbit_coupling(*ctrl, *args, x, dxdt);

            GR_precession(*ctrl, *args, x, dxdt);

            GW_radiation(*ctrl, *args, x, dxdt);
        }

        Controler const *ctrl;
        SecularConst const *args;
    };
} // namespace secular
#endif
