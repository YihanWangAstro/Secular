#ifndef DESITTER_H
#define DESITTER_H

#include "tools.h"

namespace secular {

    template<typename Container>
    struct spin_num{
        static constexpr size_t size{(std::tuple_size<Container>::value - 12)/3 };
    };

    template<bool DA, typename Args, typename Container, size_t var_idx>
    inline auto LL_coupling(Args const& args, Container const& var){
      /*---------------------------------------------------------------------------*\
          mapping alias
      \*---------------------------------------------------------------------------*/
      constexpr size_t v_offset = 3*var_idx;

      const auto[Sx, Sy, Sz] = std::tie(var[v_offset], var[v_offset + 1], var[v_offset + 2]);

      if constexpr (DA) {
          /*---------------------------------------------------------------------------*\
              mapping alias
          \*---------------------------------------------------------------------------*/
          constexpr size_t L_offset = 6;

          constexpr size_t e_offset = 9;

          const auto[Lx, Ly, Lz] = std::tie(var[L_offset], var[L_offset + 1], var[L_offset + 2]);

          const auto[ex, ey, ez] = std::tie(var[e_offset], var[e_offset + 1], var[e_offset + 2]);
          /*---------------------------------------------------------------------------*\
              orbital parameters calculation
          \*---------------------------------------------------------------------------*/
          double const a_eff = calc_a_eff(args.a_coef[1], Lx, Ly, Lz, ex, ey, ez);

          double const r3_a_eff = 1/(a_eff * a_eff * a_eff);
          /*---------------------------------------------------------------------------*\
              combinations
          \*---------------------------------------------------------------------------*/
          return cross_with_coef(args.LL * r3_a_eff, Lx, Ly, Lz, Sx, Sy, Sz);
      } else {
          const auto[rx, ry, rz] = std::tie(var[6], var[7], var[8]);

          const auto[vx, vy, vz] = std::tie(var[9], var[10], var[11]);

          auto [Lx, Ly, Lz] = cross_with_coef(args.mu[1], rx, ry, rz, vx, vy, vz);
          /*---------------------------------------------------------------------------*\
              orbital parameters calculation
          \*---------------------------------------------------------------------------*/
          double const a_eff = norm(rx, ry, rz);

          double const r3_a_eff = 1/(a_eff * a_eff * a_eff);

          return cross_with_coef(args.LL * r3_a_eff, Lx, Ly, Lz, Sx, Sy, Sz);
      }
    }

    template<bool DA, typename Args, typename Container, int S_idx, int L_idx>
    inline auto SL_coupling(Args const& args, Container const& var){
      /*---------------------------------------------------------------------------*\
          mapping alias
      \*---------------------------------------------------------------------------*/
      constexpr size_t S_offset = 12 + S_idx * 3;

      const auto[Sx, Sy, Sz] = std::tie(var[S_offset], var[S_offset + 1], var[S_offset + 2]);

      if constexpr(L_idx == 1 && !DA) {
          const auto[rx, ry, rz] = std::tie(var[6], var[7], var[8]);

          const auto[vx, vy, vz] = std::tie(var[9], var[10], var[11]);

          auto [Lx, Ly, Lz] = cross_with_coef(args.mu[1], rx, ry, rz, vx, vy, vz);
          /*---------------------------------------------------------------------------*\
              orbital parameters calculation
          \*---------------------------------------------------------------------------*/
          double const a_eff = norm(rx, ry, rz);

          double const r3_a_eff = 1/(a_eff * a_eff * a_eff);
          /*---------------------------------------------------------------------------*\
              combinations
          \*---------------------------------------------------------------------------*/
          return cross_with_coef(args.SL[S_idx][L_idx] * r3_a_eff, Lx, Ly, Lz, Sx, Sy, Sz);
      } else {
          /*---------------------------------------------------------------------------*\
              mapping alias
          \*---------------------------------------------------------------------------*/
          constexpr size_t L_offset = L_idx * 6;

          constexpr size_t e_offset = 3 + L_offset;

          const auto[Lx, Ly, Lz] = std::tie(var[L_offset], var[L_offset + 1], var[L_offset + 2]);

          const auto[ex, ey, ez] = std::tie(var[e_offset], var[e_offset + 1], var[e_offset + 2]);
          /*---------------------------------------------------------------------------*\
              orbital parameters calculation
          \*---------------------------------------------------------------------------*/
          double const a_eff = calc_a_eff(args.a_coef[L_idx], Lx, Ly, Lz, ex, ey, ez);

          double const r3_a_eff = 1/(a_eff * a_eff * a_eff);
          /*---------------------------------------------------------------------------*\
              combinations
          \*---------------------------------------------------------------------------*/
          return cross_with_coef(args.SL[S_idx][L_idx] * r3_a_eff, Lx, Ly, Lz, Sx, Sy, Sz);
      }

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
          double const r2 = norm2(rx, ry, rz);

          double const a_eff = sqrt(r2);

          double const r5_a_eff = 1/(a_eff * a_eff * a_eff * a_eff * a_eff);

          double const acc_coef = args.SL[S_idx][L_idx]*r5_a_eff;

          auto [csvx, csvy, csvz] = cross(Sx, Sy, Sz, vx, vy, vz);

          auto [csrx, csry, csrz] = cross(Sx, Sy, Sz, rx, ry, rz);

          double dvr = dot(vx, vy, vz, rx, ry, rz);

          double tri_dot = dot(Sx, Sy, Sz, crvx, crvy, crvz);

          double acc_x = acc_coef *( 3 * tri_dot * rx + 2 * r2*csvx - 3 * dvr * csrx);

          double acc_y = acc_coef *( 3 * tri_dot * ry + 2 * r2*csvy - 3 * dvr * csry);

          double acc_z = acc_coef *( 3 * tri_dot * rz + 2 * r2*csvz - 3 * dvr * csrz);

          return std::make_tuple(acc_x, acc_y, acc_z);
      } else {
          constexpr size_t L_offset = L_idx * 6;

          constexpr size_t e_offset = 3 + L_offset;

          const auto[Lx, Ly, Lz] = std::tie(var[L_offset], var[L_offset + 1], var[L_offset + 2]);

          const auto[ex, ey, ez] = std::tie(var[e_offset], var[e_offset + 1], var[e_offset + 2]);
          /*---------------------------------------------------------------------------*\
              orbital parameters calculation
          \*---------------------------------------------------------------------------*/
          auto [nex, ney, nez] = deSitter_e_vec(Sx, Sy, Sz, Lx, Ly, Lz);

          double const a_eff = calc_a_eff(args.a_coef[L_idx], Lx, Ly, Lz, ex, ey, ez);

          double const r3_a_eff = 1/(a_eff * a_eff * a_eff);
          /*---------------------------------------------------------------------------*\
              combinations
          \*---------------------------------------------------------------------------*/
          return cross_with_coef(args.SL[S_idx][L_idx] * r3_a_eff, nex, ney, nez, ex, ey, ez);
      }
    }

    template<bool DA, bool LL, bool SL, typename Args, typename Container>
    inline void deSitter_precession(Args const &args, Container const &var, Container &ddt, double t) {
        constexpr int Lin_idx = 0;
        constexpr int Lout_idx = 1;
        constexpr int S1_idx = 0;
        constexpr int S2_idx = 1;
        constexpr int S3_idx = 2;

        if constexpr(spin_num<Container>::size >= 1) {
            auto [dS1x, dS1y, dS1z] = std::tie(ddt[12], ddt[13], ddt[14]);

            std::tie(dS1x, dS1y, dS1z) =  SL_coupling<DA, Args, Container, S1_idx, Lin_idx>(args, var);
            if constexpr(SL) {
                auto [dx, dy, dz] = SL_coupling<DA, Args, Container, S1_idx, Lout_idx>(args, var);
                dS1x += dx, dS1y += dy, dS1z += dz;
            }
        }

        if constexpr(spin_num<Container>::size >= 2) {
            auto [dS2x, dS2y, dS2z] = std::tie(ddt[15], ddt[16], ddt[17]);

            std::tie(dS2x, dS2y, dS2z) = SL_coupling<DA, Args, Container, S2_idx, Lin_idx>(args, var);
            if constexpr(SL) {
                auto [dx, dy, dz] = SL_coupling<DA, Args, Container, S2_idx, Lout_idx>(args, var);
                dS2x += dx, dS2y += dy, dS2z += dz;
            }
        }

        if constexpr(spin_num<Container>::size == 3) {
            auto [dS3x, dS3y, dS3z] = std::tie(ddt[18], ddt[19], ddt[20]);

            std::tie(dS3x, dS3y, dS3z) = SL_coupling<DA, Args, Container, S3_idx, Lout_idx>(args, var);

            if constexpr(DA){
                auto [dL2x, dL2y, dL2z] = std::tie(ddt[6], ddt[7], ddt[8]);

                auto [de2x, de2y, de2z] = std::tie(ddt[9], ddt[10], ddt[11]);

                dL2x -= dS3x, dL2y -= dS3y, dL2z -= dS3z;

                auto [dex, dey, dez] = SL_coupling_bc<DA, Args, Container, S3_idx, Lout_idx>(args, var);

                de2x += dex, de2y += dey, de2z += dez;
            } else{
                auto[dvx, dvy, dvz] = std::tie(ddt[9], ddt[10], ddt[11]);

                auto [ax, ay, az] = SL_coupling_bc<DA, Args, Container, S3_idx, Lout_idx>(args, var);

                dvx += ax, dvy += ay, dvz += az;
            }
        }

        if constexpr (LL) {
            auto [dL1x, dL1y, dL1z] = std::tie(ddt[0], ddt[1], ddt[2]);

            auto [de1x, de1y, de1z] = std::tie(ddt[3], ddt[4], ddt[5]);

            auto [dLx, dLy, dLz] = LL_coupling<DA, Args, Container, 0>(args, var);//evolve L1

            auto [dex, dey, dez] = LL_coupling<DA, Args, Container, 1>(args, var);//evolve e1

            dL1x += dLx, dL1y += dLy, dL1z += dLz;

            de1x += dex, de1y += dey, de1z += dez;
        }
    }

} // namespace secular
#endif
