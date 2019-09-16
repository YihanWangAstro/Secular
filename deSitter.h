#ifndef DESITTER_H
#define DESITTER_H

#include "tools.h"

namespace secular {
    
    class SLConst {
    public:
        SLConst(double m1, double m2, double m3) {
            DS_[0] = deSitter_coef(m1, m2);
            DS_[1] = deSitter_coef(m1 + m2, m3);
            
            
            DS_[2] = deSitter_coef(m2, m1);
            DS_[3] = deSitter_coef(m1 + m2, m3);
            LT_ = 0.5 * consts::G / (consts::C * consts::C); 
            
            
            DS_[4] = deSitter_coef(m3, m1 + m2);
            DS_[5] = DS_[4];
            
            LL_ = deSitter_coef(m1 + m2, m3);
        }

        SLConst() = default;

        READ_GETTER(double, LL, LL_);

        READ_GETTER(double, S1L1, DS_[0]);

        READ_GETTER(double, S1L2, DS_[1]);

        READ_GETTER(double, S2L1, DS_[2]);

        READ_GETTER(double, S2L2, DS_[3]);

        READ_GETTER(double, S3L1, DS_[4]);

        READ_GETTER(double, S3L2, DS_[5]);

        READ_GETTER(double, S1S2, LT_);

        READ_GETTER(double, S1S3, LT_);

        READ_GETTER(double, S2S3, LT_);

    private:
        double DS_[6];
        double LT_;
        double LL_;

        inline double deSitter_coef(double m_self, double m_other) {
            return 0.5 * consts::G / (consts::C * consts::C) * (4 + 3 * m_other / m_self);
        };
    };

    enum class deS {
        off, on, bc, all
    };

    size_t to_index(deS x){
        if(x == deS::off)
            return 0;
        else if(x == deS::on)
            return 1;
        else if(x == deS::bc)
            return 2;
        else if(x == deS::all)
            return 3;
        else
            return 0;
    }

    deS str_to_spin_orbit_enum(std::string const &key)
    {
        if (case_insens_equals(key, "off") || case_insens_equals(key, "no") ){
            return deS::off;
        } else if (case_insens_equals(key, "on")) {
            return deS::on;
        } else if (case_insens_equals(key, "backreaction")) {
            return deS::bc;
        } else if (case_insens_equals(key, "both") || case_insens_equals(key, "on+backreaction")) {
            return deS::all;
        } else {
            throw ReturnFlag::input_err;
        }
    }

    template<typename Ctrl>
    bool is_Lin_coupled(Ctrl const& ctrl) {
        return ctrl.LL != deS::off || ctrl.Sin_Lin != deS::off || ctrl.Sout_Lin != deS::off || ctrl.Sin_Sin != deS::off; 
    };

    template<typename Ctrl>
    bool is_Lout_coupled(Ctrl const& ctrl) {
        return ctrl.LL != deS::off || ctrl.Sin_Lout != deS::off || ctrl.Sout_Lout != deS::off || ctrl.Sin_Sout != deS::off;
    };

    template<typename Ctrl, typename Args, typename Container>
    class deSitter_arg {
    public:
        deSitter_arg(Ctrl const& ctrl, Args const &args, Container const &var) {
            bool const  Lin_coupled{is_Lin_coupled(ctrl)};

            bool const Lout_coupled{is_Lout_coupled(ctrl)}; 

            if (Lin_coupled == true) {
                a_in_eff_ = calc_a_eff(args.a_in_coef(), var.L1x(), var.L1y(), var.L1z(), var.e1x(), var.e1y(), var.e1z());

                a_in_eff3_ = a_in_eff_ * a_in_eff_ * a_in_eff_;
            }

            if (Lout_coupled == true) {
                if (ctrl.ave_method == LK_method::DA) {
                    L2x_ = var.L2x(), L2y_ = var.L2y(), L2z_ = var.L2z();

                    a_out_eff_ = calc_a_eff(args.a_out_coef(), var.L2x(), var.L2y(), var.L2z(), var.e2x(), var.e2y(), var.e2z());

                    a_out_eff3_ = a_out_eff_ * a_out_eff_ * a_out_eff_;
                } else if(ctrl.ave_method == LK_method::SA){
                    std::tie(L2x_, L2y_, L2z_) = cross_with_coef(args.mu_out(), var.rx(), var.ry(), var.rz(), var.vx(), var.vy(), var.vz());

                    a_out_eff_ = norm(var.rx(), var.ry(), var.rz());

                    a_out_eff3_ = a_out_eff_ * a_out_eff_ * a_out_eff_;
                }
            }

            if (ctrl.Sin_Lin != deS::off)
                Omega_[0] = args.S1L1() / a_in_eff3_;
            if (ctrl.Sin_Lout != deS::off)
                Omega_[1] = args.S1L2() / a_out_eff3_;
            
            if (ctrl.Sin_Lin != deS::off)
                Omega_[2] = args.S2L1() / a_in_eff3_;
            if (ctrl.Sin_Lout != deS::off)
                Omega_[3] = args.S2L2() / a_out_eff3_;
            if (ctrl.Sin_Sin != deS::off)
                Omega_[6] = args.S1S2() / a_in_eff3_;
            
            if (ctrl.Sout_Lin != deS::off)
                Omega_[4] = args.S3L1() / a_in_eff3_;
            if (ctrl.Sout_Lout != deS::off)
                Omega_[5] = args.S3L2() / a_out_eff3_;
            if (ctrl.Sin_Sout != deS::off)
                Omega_[7] = args.S1S3() / a_out_eff3_;
            
            if (ctrl.LL != deS::off) {
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

        READ_GETTER(double, S1S2_Omega, Omega_[6]);

        READ_GETTER(double, S1S3_Omega, Omega_[7]);

        READ_GETTER(double, S2S3_Omega, Omega_[7]);

        READ_GETTER(double, L2x, L2x_);

        READ_GETTER(double, L2y, L2y_);

        READ_GETTER(double, L2z, L2z_);
    private:
        double L2x_;
        double L2y_;
        double L2z_;
        double Omega_[8];
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

    template<typename Container>
    auto SA_back_reaction(double Omega, double Sx, double Sy, double Sz, Container const &var) {
        auto[crvx, crvy, crvz] = cross(var.rx(), var.ry(), var.rz(), var.vx(), var.vy(), var.vz());

        double r2 = norm2(var.rx(), var.ry(), var.rz());

        double const acc_coef = Omega / r2;

        auto[csvx, csvy, csvz] = cross(Sx, Sy, Sz, var.vx(), var.vy(), var.vz());

        auto[csrx, csry, csrz] = cross(Sx, Sy, Sz, var.rx(), var.ry(), var.rz());

        double dvr = dot(var.rx(), var.ry(), var.rz(), var.vx(), var.vy(), var.vz());

        double tri_dot = dot(Sx, Sy, Sz, crvx, crvy, crvz);

        double acc_x = acc_coef * (3 * tri_dot * var.rx() + 2 * r2 * csvx - 3 * dvr * csrx);

        double acc_y = acc_coef * (3 * tri_dot * var.ry() + 2 * r2 * csvy - 3 * dvr * csry);

        double acc_z = acc_coef * (3 * tri_dot * var.rz() + 2 * r2 * csvz - 3 * dvr * csrz);

        return std::make_tuple(acc_x, acc_y, acc_z);
    }


#define DESITTER_IN(C, OMEGA, S)                                                                                                             \
    if(C != deS::off){                                                                                                                       \
        auto [dx, dy, dz] = cross_with_coef(OMEGA, var.L1x(), var.L1y(), var.L1z(), var.S##x(), var.S##y(), var.S##z());                     \
        if (C == deS::on || C == deS::all){                                                                                                  \
            dvar.add_##S(dx, dy, dz);                                                                                                        \
        }                                                                                                                                    \
        if (C == deS::bc || C == deS::all) {                                                                                                 \
            dvar.sub_L1(dx, dy, dz);                                                                                                         \
            auto [nex, ney, nez] = deSitter_e_vec(var.S##x(), var.S##y(), var.S##z(), var.L1x(), var.L1y(), var.L1z());                      \
            dvar.add_e1(cross_with_coef(OMEGA, nex, ney, nez, var.e1x(), var.e1y(), var.e1z()));                                             \
        }                                                                                                                                    \
    }
    


#define DESITTER_OUT(C, OMEGA, S)                                                                                                            \
    if(C != deS::off){                                                                                                                       \
        auto [dx, dy, dz] = cross_with_coef(OMEGA, d.L2x(), d.L2y(), d.L2z(), var.S##x(), var.S##y(), var.S##z());                           \
        if (C == deS::on || C == deS::all){                                                                                                  \
            dvar.add_##S(dx, dy, dz);                                                                                                        \
        }                                                                                                                                    \
        if (C == deS::bc || C == deS::all) {                                                                                                 \
            if (ctrl.ave_method == LK_method::DA){                                                                                            \
                dvar.sub_L2(dx, dy, dz);                                                                                                     \
                auto [nex, ney, nez] = deSitter_e_vec(var.S##x(), var.S##y(), var.S##z(), var.L2x(), var.L2y(), var.L2z());                  \
                dvar.add_e2(cross_with_coef(OMEGA, nex, ney, nez, var.e2x(), var.e2y(), var.e2z()));                                         \
            } else if(ctrl.ave_method == LK_method::SA) {                                                                                     \
                dvar.add_v(SA_back_reaction(OMEGA, var.S##x(), var.S##y(), var.S##z(), var));                                                \
            }                                                                                                                                \
        }                                                                                                                                    \
    }
        

#define LENS_THIRRING_IN(C, OMEGA, SI, SJ)                                                                                                   \
    if (C == deS::on || C == deS::all){                                                                                                      \
        auto [nex, ney, nez] = deSitter_e_vec(var.SI##x(), var.SI##y(), var.SI##z(), var.L1x(), var.L1y(), var.L1z());                       \
        auto [dx, dy, dz] = cross_with_coef(OMEGA, nex, ney, nez, var.SJ##x(), var.SJ##y(), var.SJ##z());                                    \
        dvar.add_##SJ(dx, dy, dz);                                                                                                           \
    }                                                                                                                                        \

#define LENS_THIRRING_OUT(C, OMEGA, SI, SJ)                                                                                                  \
    if (C == deS::on || C == deS::all){                                                                                                      \
        auto [nex, ney, nez] = deSitter_e_vec(var.SI##x(), var.SI##y(), var.SI##z(), d.L2x(), d.L2y(), d.L2z());                             \
        auto [dx, dy, dz] = cross_with_coef(OMEGA, nex, ney, nez, var.SJ##x(), var.SJ##y(), var.SJ##z());                                    \
        dvar.add_##SJ(dx, dy, dz);                                                                                                           \
    }                                                                                                                                        \


    template<typename Control, typename Args, typename Container>
    void spin_orbit_coupling(Control const& ctrl, Args const &args, Container const &var, Container &dvar) {
        using deArgs =  deSitter_arg<Control, Args, Container>;
        deArgs d{ctrl, args, var};//calculate the Omega and L2(Single average case)
        
        DESITTER_IN(ctrl.Sin_Lin, d.S1L1_Omega(), S1);
        DESITTER_OUT(ctrl.Sin_Lout, d.S1L2_Omega(), S1);
        
        DESITTER_IN(ctrl.Sin_Lin, d.S2L1_Omega(), S2);
        DESITTER_OUT(ctrl.Sin_Lout, d.S2L2_Omega(), S2);

        LENS_THIRRING_IN(ctrl.Sin_Sin, d.S1S2_Omega(), S1, S2);
        LENS_THIRRING_IN(ctrl.Sin_Sin, d.S1S2_Omega(), S2, S1);
        
        LENS_THIRRING_OUT(ctrl.Sout_Lin, d.S3L1_Omega(), S3, L1);
        LENS_THIRRING_OUT(ctrl.Sout_Lin, d.S3L1_Omega(), S3, e1);//if the back rection on, L1 E1 will back react to Lout twice !!need to fix it
        LENS_THIRRING_OUT(ctrl.Sout_Lin, d.S3L1_Omega(), L1, S1);
       
        DESITTER_OUT(ctrl.Sout_Lout, d.S3L2_Omega(), S3);

        LENS_THIRRING_OUT(ctrl.Sin_Sout, d.S1S3_Omega(), S3, S1);
        LENS_THIRRING_OUT(ctrl.Sin_Sout, d.S1S3_Omega(), S1, S3);
        LENS_THIRRING_OUT(ctrl.Sin_Sout, d.S2S3_Omega(), S3, S2);
        LENS_THIRRING_OUT(ctrl.Sin_Sout, d.S2S3_Omega(), S2, S3);

        DESITTER_OUT(ctrl.LL, d.LL(), L1);
        DESITTER_OUT(ctrl.LL, d.LL(), e1);
    }

} // namespace secular
#endif
