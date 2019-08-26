#ifndef DESITTER_H
#define DESITTER_H

#include "tools.h"

namespace secular {

    template<typename Container>
    struct spin_num{
        static constexpr size_t size{(std::tuple_size<Container>::value - 12)/3 };
    };

    template<size_t sp_num>
    class SLConst{
        static_assert(spin_num<=3, "The max spin number is 3!");
    public:
        SLConst(double m1, double m2, double m3) {
            if constexpr(sp_num > 0) {
                SL_[0] = deSitter_coef(m1, m2);
                SL_[1] = deSitter_coef(m1+m2, m3);
            }
            if constexpr(sp_num > 1) {
                SL_[2] = deSitter_coef(m2, m1);
                SL_[3] = deSitter_coef(m1+m2, m3);
            }
            if constexpr(sp_num == 3) {
                SL_[4] = deSitter_coef(m3, m1+m2);
                SL_[5] = SL_[4];
            }
            LL_ = deSitter_coef(m1+m2, m3);
        }

        SLConst() = default;

        READ_GETTER(double, LL, LL_);
        OPT_READ_GETTER(sp_num>0, double, S1L1, SL_[0]);
        OPT_READ_GETTER(sp_num>0, double, S1L2, SL_[1]);
        OPT_READ_GETTER(sp_num>1, double, S2L1, SL_[2]);
        OPT_READ_GETTER(sp_num>1, double, S2L2, SL_[3]);
        OPT_READ_GETTER(sp_num>2, double, S3L1, SL_[4]);
        OPT_READ_GETTER(sp_num>2, double, S3L2, SL_[5]);
    private:
        double SL_[sp_num*2];
        double LL_;

        inline double deSitter_coef(double m_self, double m_other) {
            return 0.5*consts::G/(consts::C * consts::C) * (4 + 3*m_other/m_self);
        };
    };



    template<bool DA, bool LL, bool SL, typename Args, typename Container>
    struct deSitter_arg{

        deSitter_arg(Args const&args, Container const& var) {
            Lx[0] = var.L1x(), Ly[0] = var.L1y(), Lz[0] = var.L1z();

            a_eff[0] = calc_a_eff(args.a_in_coef(), var.L1x(), var.L1y(), var.L1z(), var.e1x(), var.e1y(), var.e1z());

            a_eff3[0] = a_eff[0]* a_eff[0]*a_eff[0];

            if constexpr(LL || SL || spin_num<Container>::size == 3){
                if constexpr(DA){
                    Lx[1] = var.L2x(), Ly[1] = var.L2y(), Lz[1] = var.L2z();

                    a_eff[1] = calc_a_eff(args.a_out_coef(), var.L2x(), var.L2y(), var.L2z(), var.e2x(), var.e2y(), var.e2z());

                    a_eff3[1] = a_eff[1]* a_eff[1]*a_eff[1];
                } else {
                    std::tie(Lx[1], Ly[1], Lz[1]) = cross_with_coef(args.mu_in(),var.rx(), var.ry(), var.rz(), var.vx(), var.vy(), var.vz());

                    a_eff[1] = norm(var.rx(), var.ry(), var.rz());

                    a_eff3[1] = a_eff[1]* a_eff[1]*a_eff[1];
                }
            }

            if constexpr(spin_num<Container>::size >= 1){
                Omega[0][0] = args.SL(0,0) / a_eff3[0];
                if constexpr(SL)
                    Omega[0][1] = args.SL(0,1) / a_eff3[1];
            }

            if constexpr(spin_num<Container>::size >= 2){
                Omega[1][0] = args.SL(1, 0) / a_eff3[0];
                if constexpr(SL)
                    Omega[1][1] = args.SL(1,1) / a_eff3[1];
            }

            if constexpr(spin_num<Container>::size == 3){
                Omega[2][0] = args.SL(2,1) / a_eff3[1];
                Omega[2][1] = Omega[2][0];

            }
            if constexpr(LL) {
                LL_Omega = args.LL() / a_eff3[1];
            }
        }

        double Lx[2];
        double Ly[2];
        double Lz[2];
        double Omega[3][2];
        double LL_Omega{0};
        double a_eff[2];
        double a_eff3[2];
    };

    inline auto deSitter_e_vec(double S1x, double S1y, double S1z, double Lx, double Ly, double Lz) {
        double dot_part = 3*dot(Lx, Ly, Lz, S1x, S1y, S1z)/norm2(Lx, Ly, Lz);
        return std::make_tuple(S1x - dot_part*Lx, S1y - dot_part*Ly, S1z - dot_part*Lz);
    }

    template<bool DA, typename Args, typename Container, int S_idx, int L_idx>
    inline auto SL_coupling_bc(Args const& args, Container const& var){
      /*---------------------------------------------------------------------------*\
          mapping alias
      \*---------------------------------------------------------------------------*/
      constexpr size_t S_offset = 12 + S_idx * 3;

      const auto[Sx, Sy, Sz] = std::tie(var[S_offset], var[S_offset + 1], var[S_offset + 2]);

      if constexpr(L_idx == 1 && !DA) {
          auto [crvx, crvy, crvz] = cross(var.rx(), var.ry(), var.rz(), var.vx(), var.vy(), var.vz());
          /*---------------------------------------------------------------------------*\
              orbital parameters calculation
          \*---------------------------------------------------------------------------*/
          double r2 = norm2(var.rx(), var.ry(), var.rz());

          double const acc_coef = args.Omega[S_idx][L_idx] / r2;

          auto [csvx, csvy, csvz] = cross(Sx, Sy, Sz, var.vx(), var.vy(), var.vz());

          auto [csrx, csry, csrz] = cross(Sx, Sy, Sz, var.rx(), var.ry(), var.rz());

          double dvr = dot(var.rx(), var.ry(), var.rz(), var.vx(), var.vy(), var.vz());

          double tri_dot = dot(Sx, Sy, Sz, crvx, crvy, crvz);

          double acc_x = acc_coef *( 3 * tri_dot * var.rx() + 2 * r2*csvx - 3 * dvr * csrx);

          double acc_y = acc_coef *( 3 * tri_dot * var.ry() + 2 * r2*csvy - 3 * dvr * csry);

          double acc_z = acc_coef *( 3 * tri_dot * var.rz() + 2 * r2*csvz - 3 * dvr * csrz);

          return std::make_tuple(acc_x, acc_y, acc_z);
      } else {
          constexpr size_t e_offset = 3 + L_idx * 6;

          const auto[ex, ey, ez] = std::tie(var[e_offset], var[e_offset + 1], var[e_offset + 2]);

          auto [nex, ney, nez] = deSitter_e_vec(Sx, Sy, Sz, args.Lx[L_idx], args.Ly[L_idx], args.Lz[L_idx]);

          return cross_with_coef(args.Omega[S_idx][L_idx], nex, ney, nez, ex, ey, ez);
      }
    }

    template<bool DA, bool LL, bool SL, typename Args, typename Container>
    inline void deSitter_precession(Args const &args, Container const &var, Container &dvar, double t) {
        constexpr int Lin_idx = 0;
        constexpr int Lout_idx = 1;
        constexpr int S1_idx = 0;
        constexpr int S2_idx = 1;
        constexpr int S3_idx = 2;
        using deArgs =  deSitter_arg<DA, LL, SL, Args, Container>;
        deArgs d{args, var};

        if constexpr(spin_num<Container>::size >= 1) {
            dvar.set_S1(cross_with_coef(d.Omega(S1_idx, Lin_idx), var.L1x(), var.L1y(), var.L1z(), var.S1x(), var.S1y(), var.S1z()));
            if constexpr(SL) {
                dvar.add_S1(cross_with_coef(d.Omega(S1_idx, Lout_idx), d.L2x(), d.L2y(), d.L2z(), var.S1x(), var.S1y(), var.S1z()));
            }
        }

        if constexpr(spin_num<Container>::size >= 2) {
            dvar.set_S2(cross_with_coef(d.Omega(S2_idx, Lin_idx), var.L1x(), var.L1y(), var.L1z(), var.S2x(), var.S2y(), var.S2z()));
            if constexpr(SL) {
                dvar.add_S2(cross_with_coef(d.Omega(S2_idx, Lout_idx), d.L2x(), d.L2y(), d.L2z(), var.S2x(), var.S2y(), var.S2z()));
            }
        }

        if constexpr(spin_num<Container>::size == 3) {
            dvar.set_S3(cross_with_coef(d.Omega(S3_idx, Lout_idx), d.L2x(), d.L2y(), d.L2z(), var.S3x(), var.S3y(), var.S3z()));
            if constexpr(DA){
                dvar.sub_L2(dvar.S3x(), dvar.S3y(), dvar.S3z())
                dvar.add_e2(SL_coupling_bc<DA, deArgs, Container, S3_idx, Lout_idx>(deS_args, var));
            } else{
                dvar.add_v(SL_coupling_bc<DA, deArgs, Container, S3_idx, Lout_idx>(deS_args, var));
            }
        }

        if constexpr (LL) {
            dvar.add_L1(cross_with_coef(d.LL(), d.L2x(), d.L2y(), d.L2z(), var.L1x(), var.L1y(), var.L1z()));//evolve L1
            dvar.add_e1(cross_with_coef(d.LL(), d.L2x(), d.L2y(), d.L2z(), var.e1x(), var.e1y(), var.e1z()));//evolve e1
        }
    }

} // namespace secular
#endif
