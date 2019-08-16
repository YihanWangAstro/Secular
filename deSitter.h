#ifndef DESITTER_H
#define DESITTER_H

#include "tools.h"

namespace secular {

    template<typename Container>
    constexpr size_t spin_num(Container const &d) {
        if constexpr (d.size() == 12) {
            return 0;
        } else if constexpr (d.size() == 15) {
            return 1;
        } else if constexpr (d.size() == 18) {
            return 2;
        } else if constexpr (d.size() == 21) {
            return 3;
        } else {
            std::cout << "Wrong Spin Num!\n";
            exit(0);
        }
    }

    template<bool DA, typename Args, typename Container, int S_idx, int L_idx>
    inline auto SL_coupling(Args const& args, Container const& var){
      /*---------------------------------------------------------------------------*\
          mapping alias
      \*---------------------------------------------------------------------------*/
      constexpr size_t S_offset = 12 + s_idx * 3;

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

    template<bool DA, typename Args, typename Container, int S_idx, int L_idx>
    inline auto SL_coupling_bc(Args const& args, Container const& var){
      /*---------------------------------------------------------------------------*\
          mapping alias
      \*---------------------------------------------------------------------------*/
      constexpr size_t S_offset = 12 + s_idx * 3;

      const auto[Sx, Sy, Sz] = std::tie(var[S_offset], var[S_offset + 1], var[S_offset + 2]);

      if constexpr(L_idx == 1 && !DA) {
          /*---------------------------------------------------------------------------*\
              mapping alias
          \*---------------------------------------------------------------------------*/
          const auto[rx, ry, rz] = std::tie(var[6], var[7], var[8]);

          const auto[vx, vy, vz] = std::tie(var[9], var[10], var[11]);

          auto [Lx, Ly, Lz] = cross_with_coef(args.mu[1], rx, ry, rz, vx, vy, vz);
          /*---------------------------------------------------------------------------*\
              orbital parameters calculation
          \*---------------------------------------------------------------------------*/
          double const a_eff = norm(rx, ry, rz);

          double const r3_a_eff = 1/(a_eff * a_eff * a_eff);
          ////............
          return ;
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

        if constexpr (LL) {
            SL_coupling<DA, Container, -4, Lout_idx>(args, var);//evolve L1
            SL_coupling<DA, Container, -3, Lout_idx>(args, var);//evolve e1
        }

        if constexpr(spin_num(var) >= 1) {
            std::tie(var[S_offset], var[S_offset + 1], var[S_offset + 2]) =  SL_coupling<DA, Container, S1_idx, Lin_idx>(args, var);
            if constexpr(SL) {
                SL_coupling<DA, Container, S1_idx, Lout_idx>(args, var);
            }
        }

        if constexpr(spin_num(var) >= 2) {
            SL_coupling<DA, Container, S2_idx, Lin_idx>(args, var);
            if constexpr(SL) {
                SL_coupling<DA, Container, S2_idx, Lout_idx>(args, var);
            }
        }

        if constexpr(spin_num(var) == 3) {
            SL_coupling<DA, Container, S3_idx, L_out_idx>(args, var);
        }
    }

} // namespace secular
#endif
