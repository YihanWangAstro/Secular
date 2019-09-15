#ifndef DESITTER_H
#define DESITTER_H

#include "tools.h"

namespace secular {
    template<size_t sp_num>
    class SLConst {
        static_assert(sp_num <= 3, "The max spin number is 3!");
    public:
        SLConst(double m1, double m2, double m3) {
            if constexpr(sp_num > 0) {
                SL_[0] = deSitter_coef(m1, m2);
                SL_[1] = deSitter_coef(m1 + m2, m3);
            }
            if constexpr(sp_num > 1) {
                SL_[2] = deSitter_coef(m2, m1);
                SL_[3] = deSitter_coef(m1 + m2, m3);
            }
            if constexpr(sp_num == 3) {
                SL_[4] = deSitter_coef(m3, m1 + m2);
                SL_[5] = SL_[4];
            }
            LL_ = deSitter_coef(m1 + m2, m3);
        }

        SLConst() = default;

        READ_GETTER(double, LL, LL_);

        OPT_READ_GETTER(sp_num > 0, double, S1L1, SL_[0]);

        OPT_READ_GETTER(sp_num > 0, double, S1L2, SL_[1]);

        OPT_READ_GETTER(sp_num > 1, double, S2L1, SL_[2]);

        OPT_READ_GETTER(sp_num > 1, double, S2L2, SL_[3]);

        OPT_READ_GETTER(sp_num > 2, double, S3L1, SL_[4]);

        OPT_READ_GETTER(sp_num > 2, double, S3L2, SL_[5]);
    private:
        double SL_[sp_num * 2+1];
        double LL_;

        inline double deSitter_coef(double m_self, double m_other) {
            return 0.5 * consts::G / (consts::C * consts::C) * (4 + 3 * m_other / m_self);
        };
    };

    enum class deS {
        off, on, bc, all
    };

    template<typename Stat>
    struct is_Lin_coupled {
        constexpr static bool value{Stat::LL != deS::off || Stat::Sin_Lin != deS::off || Stat::Sout_Lin != deS::off};
    };

    template<typename Stat>
    struct is_Lout_coupled {
        constexpr static bool value{Stat::LL != deS::off || Stat::Sin_Lout != deS::off || Stat::Sout_Lout != deS::off};
    };

    template<bool DA, typename Stat, typename Args, typename Container>
    class deSitter_arg {
        static_assert(spin_num<Container>::size <= 3, "The max spin number is 3!");
    public:
        deSitter_arg(Args const &args, Container const &var) {
            if constexpr(spin_num<Container>::size > 0 && Lin_coupled) {
                a_in_eff_ = calc_a_eff(args.a_in_coef(), var.L1x(), var.L1y(), var.L1z(), var.e1x(), var.e1y(), var.e1z());

                a_in_eff3_ = a_in_eff_ * a_in_eff_ * a_in_eff_;
            }

            if constexpr(spin_num<Container>::size > 0 && Lout_coupled) {
                if constexpr(DA) {
                    L2x_ = var.L2x(), L2y_ = var.L2y(), L2z_ = var.L2z();

                    a_out_eff_ = calc_a_eff(args.a_out_coef(), var.L2x(), var.L2y(), var.L2z(), var.e2x(), var.e2y(), var.e2z());

                    a_out_eff3_ = a_out_eff_ * a_out_eff_ * a_out_eff_;
                } else {
                    std::tie(L2x_, L2y_, L2z_) = cross_with_coef(args.mu_out(), var.rx(), var.ry(), var.rz(), var.vx(), var.vy(), var.vz());

                    a_out_eff_ = norm(var.rx(), var.ry(), var.rz());

                    a_out_eff3_ = a_out_eff_ * a_out_eff_ * a_out_eff_;
                }
            }

            if constexpr(spin_num<Container>::size > 0) {
                if constexpr(Stat::Sin_Lin != deS::off)
                    Omega_[0] = args.S1L1() / a_in_eff3_;
                if constexpr(Stat::Sin_Lout != deS::off)
                    Omega_[1] = args.S1L2() / a_out_eff3_;
            }

            if constexpr(spin_num<Container>::size > 1) {
                if constexpr(Stat::Sin_Lin != deS::off)
                    Omega_[2] = args.S2L1() / a_in_eff3_;
                if constexpr(Stat::Sin_Lout != deS::off)
                    Omega_[3] = args.S2L2() / a_out_eff3_;
            }

            if constexpr(spin_num<Container>::size > 2) {
                if constexpr(Stat::Sout_Lin != deS::off)
                    Omega_[4] = args.S3L1() / a_in_eff3_;
                if constexpr(Stat::Sout_Lout != deS::off)
                    Omega_[5] = args.S3L2() / a_out_eff3_;
            }

            if constexpr(Stat::LL != deS::off) {
                LL_ = args.LL() / a_out_eff3_;
            }
        }

        constexpr static bool Lin_coupled{is_Lin_coupled<Stat>::value};

        constexpr static bool Lout_coupled{is_Lout_coupled<Stat>::value};

        READ_GETTER(double, LL, LL_);

        OPT_READ_GETTER(spin_num<Container>::size > 0, double, S1L1_Omega, Omega_[0]);

        OPT_READ_GETTER(spin_num<Container>::size > 0, double, S1L2_Omega, Omega_[1]);

        OPT_READ_GETTER(spin_num<Container>::size > 1, double, S2L1_Omega, Omega_[2]);

        OPT_READ_GETTER(spin_num<Container>::size > 1, double, S2L2_Omega, Omega_[3]);

        OPT_READ_GETTER(spin_num<Container>::size > 2, double, S3L1_Omega, Omega_[4]);

        OPT_READ_GETTER(spin_num<Container>::size > 2, double, S3L2_Omega, Omega_[5]);


        READ_GETTER(double, L2x, L2x_);

        READ_GETTER(double, L2y, L2y_);

        READ_GETTER(double, L2z, L2z_);
    private:
        double L2x_;
        double L2y_;
        double L2z_;
        double Omega_[spin_num<Container>::size * 2 + 1];
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


#define EVOLVE_SEL_IN(C, S)                                                                                                                  \
    if (C == deS::on || C == deS::all){                                                                                                      \
        auto [dx, dy, dz] = cross_with_coef(d.S##L1_Omega(), var.L1x(), var.L1y(), var.L1z(), var.S##x(), var.S##y(), var.S##z());           \
        dvar.add_##S(dx, dy, dz);                                                                                                            \
    }                                                                                                                                        \
    if (C == deS::bc || C == deS::all) {                                                                                                     \
        dvar.sub_L1(dx, dy, dz);                                                                                                             \
        auto [nex, ney, nez] = deSitter_e_vec(var.S##x(), var.S##y(), var.S##z(), var.L1x(), var.L1y(), var.L1z());                          \
        dvar.add_e1(cross_with_coef(d.S##L1_Omega(), nex, ney, nez, var.e1x(), var.e1y(), var.e1z()));                                       \
    }                                                                                                                                        \


#define EVOLVE_SEL_OUT(C, S)                                                                                                                 \
        if (C == deS::on || C == deS::all){                                                                                                  \
            auto [dx, dy, dz] = cross_with_coef(d.S##L2_Omega(), d.L2x(), d.L2y(), d.L2z(), var.S##x(), var.S##y(), var.S##z());             \
            dvar.add_##S(dx, dy, dz);                                                                                                        \
        }                                                                                                                                    \
        if (C == deS::bc || C == deS::all) {                                                                                                 \
            if constexpr (DA){                                                                                                               \
                dvar.sub_L2(dx, dy, dz);                                                                                                     \
                auto [nex, ney, nez] = deSitter_e_vec(var.S##x(), var.S##y(), var.S##z(), d.L2x(), d.L2y(), d.L2z());                        \
                dvar.add_e2(cross_with_coef(d.S##L2_Omega(), nex, ney, nez, var.e2x(), var.e2y(), var.e2z()));                               \
            } else {                                                                                                                         \
                dvar.add_v(SA_back_reaction(d.S##L2_Omega(), var.S##x(), var.S##y(), var.S##z(), var));                                      \
            }                                                                                                                                \
        }                                                                                                                                    \

#define EVOLVE_SS(C, SI, SJ, L)                                                                                                              \
    if (C == deS::on || C == deS::all){                                                                                                      \
        auto [nex, ney, nez] = deSitter_e_vec(var.SI##x(), var.SI##y(), var.SI##z(), var.L##x(), var.L##y(), var.L##z());                    \
        auto [dx, dy, dz] = cross_with_coef(d.SI##SJ##_Omega(), nex, ney, nez, var.SJ##x(), var.SJ##y(), var.SJ##z());                       \
        dvar.add_##SJ(dx, dy, dz);                                                                                                           \
    }                                                                                                                                        \


    template<bool DA, typename Control, typename Args, typename Container>
    inline void deSitter_precession(Control const& ctrl, Args const &args, Container const &var, Container &dvar, double t) {
        using deArgs =  deSitter_arg<DA, Stat, Args, Container>;
        deArgs d{args, var};//calculate the Omega and L2(Single average case)

        if constexpr(spin_num<Container>::size >= 1) {
            EVOLVE_SEL_IN(ctrl.Sin_Lin, S1);
            EVOLVE_SEL_OUT(ctrl.Sin_Lout, S1);
        }

        if constexpr(spin_num<Container>::size >= 2) {
            EVOLVE_SEL_IN(ctrl.Sin_Lin, S2);
            EVOLVE_SEL_OUT(ctrl.Sin_Lout, S2);
        }

        if constexpr(spin_num<Container>::size == 3) {
            EVOLVE_SEL_IN(ctrl.Sout_Lin, S3);
            EVOLVE_SEL_OUT(ctrl.Sout_Lout, S3);
        }

        if (ctrl.LL == deS::on || ctrl.LL == deS::all) {
            dvar.add_L1(cross_with_coef(d.LL(), d.L2x(), d.L2y(), d.L2z(), var.L1x(), var.L1y(), var.L1z()));//evolve L1
            dvar.add_e1(cross_with_coef(d.LL(), d.L2x(), d.L2y(), d.L2z(), var.e1x(), var.e1y(), var.e1z()));//evolve e1
        }
    }

} // namespace secular
#endif
