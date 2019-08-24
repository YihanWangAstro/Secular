#ifndef DESITTER_H
#define DESITTER_H

#include "tools.h"

namespace secular {

    template<typename Container>
    struct spin_num{
        static constexpr size_t size{(std::tuple_size<Container>::value - 12)/3 };
    };

    template<bool DA, bool LL, bool SL, typename Args, typename Container>
    struct deSitter_arg{

        deSitter_arg(Args const&args, Container const& var) {
            const auto[L1x, L1y, L1z] = std::tie(var[0], var[1], var[2]);

            const auto[e1x, e1y, e1z] = std::tie(var[3], var[4], var[5]);

            Lx[0] = L1x, Ly[0] = L1y, Lz[0] = L1z;

            a_eff[0] = calc_a_eff(args.a_coef[0], L1x, L1y, L1z, e1x, e1y, e1z);

            a_eff3[0] = a_eff[0]* a_eff[0]*a_eff[0];

            if constexpr(LL || SL || spin_num<Container>::size == 3){
                if constexpr(DA){
                    const auto[L2x, L2y, L2z] = std::tie(var[6], var[7], var[8]);

                    const auto[e2x, e2y, e2z] = std::tie(var[9], var[10], var[11]);

                    Lx[1] = L2x, Ly[1] = L2y, Lz[1] = L2z;

                    a_eff[1] = calc_a_eff(args.a_coef[1], L2x, L2y, L2z, e2x, e2y, e2z);

                    a_eff3[1] = a_eff[1]* a_eff[1]*a_eff[1];
                } else {
                    const auto[rx, ry, rz] = std::tie(var[6], var[7], var[8]);

                    const auto[vx, vy, vz] = std::tie(var[9], var[10], var[11]);

                    std::tie(Lx[1], Ly[1], Lz[1]) = cross_with_coef(args.mu[1], rx, ry, rz, vx, vy, vz);

                    a_eff[1] = norm(rx, ry, rz);

                    a_eff3[1] = a_eff[1]* a_eff[1]*a_eff[1];
                }
            }

            if constexpr(spin_num<Container>::size >= 1){
                Omega[0][0] = args.SL[0][0] / a_eff3[0];
                if constexpr(SL)
                    Omega[0][1] = args.SL[0][1] / a_eff3[1];
            }

            if constexpr(spin_num<Container>::size >= 2){
                Omega[1][0] = args.SL[1][0] / a_eff3[0];
                if constexpr(SL)
                    Omega[1][1] = args.SL[1][1] / a_eff3[1];
            }

            if constexpr(spin_num<Container>::size == 3){
                Omega[2][1] = args.SL[2][1] / a_eff3[1];
                Omega[2][0] = 0;
            }
            if constexpr(LL) {
                LL_Omega = args.LL / a_eff3[1];
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


    template<bool DA, typename Args, typename Container, size_t var_idx>
    inline auto LL_coupling(Args const& args, Container const& var){
      /*---------------------------------------------------------------------------*\
          mapping alias
      \*---------------------------------------------------------------------------*/
      constexpr size_t v_offset = 3*var_idx;

      const auto[Sx, Sy, Sz] = std::tie(var[v_offset], var[v_offset + 1], var[v_offset + 2]);

      return cross_with_coef(args.LL_Omega, args.Lx[1], args.Ly[1], args.Lz[1], Sx, Sy, Sz);
    }

    template<bool DA, typename Args, typename Container, int S_idx, int L_idx>
    inline auto SL_coupling(Args const& args, Container const& var){
      /*---------------------------------------------------------------------------*\
          mapping alias
      \*---------------------------------------------------------------------------*/
      constexpr size_t S_offset = 12 + S_idx * 3;

      const auto[Sx, Sy, Sz] = std::tie(var[S_offset], var[S_offset + 1], var[S_offset + 2]);

      return cross_with_coef(args.Omega[S_idx][L_idx], args.Lx[L_idx], args.Ly[L_idx], args.Lz[L_idx], Sx, Sy, Sz);
    }

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
          /*---------------------------------------------------------------------------*\
              mapping alias
          \*---------------------------------------------------------------------------*/
          const auto[rx, ry, rz] = std::tie(var[6], var[7], var[8]);

          const auto[vx, vy, vz] = std::tie(var[9], var[10], var[11]);

          auto [crvx, crvy, crvz] = cross(rx, ry, rz, vx, vy, vz);
          /*---------------------------------------------------------------------------*\
              orbital parameters calculation
          \*---------------------------------------------------------------------------*/
          double r2 = norm2(rx, ry, rz);

          double const acc_coef = args.Omega[S_idx][L_idx] / r2;

          auto [csvx, csvy, csvz] = cross(Sx, Sy, Sz, vx, vy, vz);

          auto [csrx, csry, csrz] = cross(Sx, Sy, Sz, rx, ry, rz);

          double dvr = dot(vx, vy, vz, rx, ry, rz);

          double tri_dot = dot(Sx, Sy, Sz, crvx, crvy, crvz);

          double acc_x = acc_coef *( 3 * tri_dot * rx + 2 * r2*csvx - 3 * dvr * csrx);

          double acc_y = acc_coef *( 3 * tri_dot * ry + 2 * r2*csvy - 3 * dvr * csry);

          double acc_z = acc_coef *( 3 * tri_dot * rz + 2 * r2*csvz - 3 * dvr * csrz);

          return std::make_tuple(acc_x, acc_y, acc_z);
      } else {
          constexpr size_t e_offset = 3 + L_idx * 6;

          const auto[ex, ey, ez] = std::tie(var[e_offset], var[e_offset + 1], var[e_offset + 2]);

          auto [nex, ney, nez] = deSitter_e_vec(Sx, Sy, Sz, args.Lx[L_idx], args.Ly[L_idx], args.Lz[L_idx]);

          return cross_with_coef(args.Omega[S_idx][L_idx], nex, ney, nez, ex, ey, ez);
      }
    }

    template<bool DA, bool LL, bool SL, typename Args, typename Container>
    inline void deSitter_precession(Args const &args, Container const &var, Container &ddt, double t) {
        constexpr int Lin_idx = 0;
        constexpr int Lout_idx = 1;
        constexpr int S1_idx = 0;
        constexpr int S2_idx = 1;
        constexpr int S3_idx = 2;
        using deArgs =  deSitter_arg<DA, LL, SL, Args, Container>;
        deArgs deS_args{args, var};

        if constexpr(spin_num<Container>::size >= 1) {
            auto [dS1x, dS1y, dS1z] = std::tie(ddt[12], ddt[13], ddt[14]);

            std::tie(dS1x, dS1y, dS1z) =  SL_coupling<DA, deArgs, Container, S1_idx, Lin_idx>(deS_args, var);
            if constexpr(SL) {
                auto [dx, dy, dz] = SL_coupling<DA, deArgs, Container, S1_idx, Lout_idx>(deS_args, var);
                dS1x += dx, dS1y += dy, dS1z += dz;
            }
        }

        if constexpr(spin_num<Container>::size >= 2) {
            auto [dS2x, dS2y, dS2z] = std::tie(ddt[15], ddt[16], ddt[17]);

            std::tie(dS2x, dS2y, dS2z) = SL_coupling<DA, deArgs, Container, S2_idx, Lin_idx>(deS_args, var);
            if constexpr(SL) {
                auto [dx, dy, dz] = SL_coupling<DA, deArgs, Container, S2_idx, Lout_idx>(deS_args, var);
                dS2x += dx, dS2y += dy, dS2z += dz;
            }
        }

        if constexpr(spin_num<Container>::size == 3) {
            auto [dS3x, dS3y, dS3z] = std::tie(ddt[18], ddt[19], ddt[20]);

            std::tie(dS3x, dS3y, dS3z) = SL_coupling<DA, deArgs, Container, S3_idx, Lout_idx>(deS_args, var);

            if constexpr(DA){
                auto [dL2x, dL2y, dL2z] = std::tie(ddt[6], ddt[7], ddt[8]);

                auto [de2x, de2y, de2z] = std::tie(ddt[9], ddt[10], ddt[11]);

                dL2x -= dS3x, dL2y -= dS3y, dL2z -= dS3z;

                auto [dex, dey, dez] = SL_coupling_bc<DA, deArgs, Container, S3_idx, Lout_idx>(deS_args, var);

                de2x += dex, de2y += dey, de2z += dez;
            } else{
                auto[dvx, dvy, dvz] = std::tie(ddt[9], ddt[10], ddt[11]);

                auto [ax, ay, az] = SL_coupling_bc<DA, deArgs, Container, S3_idx, Lout_idx>(deS_args, var);

                dvx += ax, dvy += ay, dvz += az;
            }
        }

        if constexpr (LL) {
            auto [dL1x, dL1y, dL1z] = std::tie(ddt[0], ddt[1], ddt[2]);

            auto [de1x, de1y, de1z] = std::tie(ddt[3], ddt[4], ddt[5]);

            auto [dLx, dLy, dLz] = LL_coupling<DA, deArgs, Container, 0>(deS_args, var);//evolve L1

            auto [dex, dey, dez] = LL_coupling<DA, deArgs, Container, 1>(deS_args, var);//evolve e1

            dL1x += dLx, dL1y += dLy, dL1z += dLz;

            de1x += dex, de1y += dey, de1z += dez;
        }
    }

} // namespace secular
#endif
