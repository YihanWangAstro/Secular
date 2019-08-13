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

    inline auto deSitter_e_vec(double v1x, double v1y, double v1z, double Lx, double Ly, double Lz) {
        double dot_part = 3*dot(Lx, Ly, Lz, v1x, v1y, v1z)/norm2(Lx, Ly, Lz);
        return std::make_tuple(v1x - dot_part*Lx, v1y - dot_part*Ly, v1z - dot_part*Lz);
    }

    template<bool BackReaction, typename Container, int s_idx, int L_idx>
    void SL_coupling(double Omega, Container const& var, Container& ddt) {
        /*---------------------------------------------------------------------------*\
            mapping alias
        \*---------------------------------------------------------------------------*/
        constexpr size_t L_offset = L_idx*6;

        constexpr size_t S_offset = 12 + s_idx*3;

        const auto[Lx, Ly, Lz] = std::tie(var[L_offset], var[L_offset + 1], var[L_offset + 2]);

        const auto[Sx, Sy, Sz] = std::tie(var[S_offset], var[S_offset + 1], var[S_offset + 2]);

        auto[dSx, dSy, dSz] = std::tie(ddt[S_offset], ddt[S_offset + 1], ddt[S_offset + 2]);
        /*---------------------------------------------------------------------------*\
            cross production
        \*---------------------------------------------------------------------------*/
        auto const[dSx__, dSy__, dSz__] = cross_with_coef(Omega, Lx, Ly, Lz, Sx, Sy, Sz);
        /*---------------------------------------------------------------------------*\
            combinations
        \*---------------------------------------------------------------------------*/
        dSx += dSx__, dSy += dSy__, dSz += dSz__;

        if constexpr(BackReaction) {
          constexpr size_t e_offset = 3 + L_offset;

          const auto[ex, ey, ez] = std::tie(var[e_offset], var[e_offset + 1], var[e_offset + 2]);

          auto[dLx, dLy, dLz] = std::tie(ddt[L_offset], ddt[L_offset + 1], ddt[L_offset + 2]);

          auto[dex, dey, dez] = std::tie(ddt[e_offset], ddt[e_offset + 1], ddt[e_offset + 2]);

          dLx -= dSx__, dLy -= dSy__, dLz -= dSz__;

          auto [nex, ney, nez] = deSitter_e_vec(Sx, Sy, Sz, Lx, Ly, Lz);

          auto const[dex__, dey__, dez__] = cross_with_coef(Omega, nex, ney, nez, ex, ey, ez);

          dex += dex__, dey += dey__, dez += dez__;
        }
    }

    template<bool DA, bool LL, bool SL, typename Args, typename Container>
    inline void deSitter_precession(Args const &args, Container const &var, Container &ddt, double t) {
        /*---------------------------------------------------------------------------*\
            mapping alias
        \*---------------------------------------------------------------------------*/
        const auto[L1x, L1y, L1z] = std::tie(var[0], var[1], var[2]);

        const auto[e1x, e1y, e1z] = std::tie(var[3], var[4], var[5]);

        double L2x, L2y, L2z;
        /*---------------------------------------------------------------------------*\
            orbital parameters calculation
        \*---------------------------------------------------------------------------*/
        double a_in_eff = calc_a_eff(args.a_in_coef, L1x, L1y, L1z, e1x, e1y, e1z);

        double r3_a_in_eff = 1/(a_in_eff*a_in_eff*a_in_eff);

        double r3_a_out_eff{0};

        if constexpr(LL || SL || spin_num(var) == 3) {//if out orbit is coupled, we need to calculate a_out_eff
            if constexpr(DA) {//double averaged
                L2x = var[6], L2y = var[7], L2z = var[8];

                const auto[e2x, e2y, e2z] = std::tie(var[9], var[10], var[11]);

                double a_out_eff = calc_a_eff(args.a_out_coef, L2x, L2y, L2z, e2x, e2y, e2z);

                r3_a_out_eff = 1/(a_out_eff*a_out_eff*a_out_eff);
            } else {//single averaged
                const auto[rx, ry, rz] = std::tie(var[6], var[7], var[8]);

                const auto[vx, vy, vz] = std::tie(var[9], var[10], var[11]);

                std::tie(L2x, L2y, L2z) = cross_with_coef(args.mu2, rx, ry, rz, vx, vy, vz);

                double a_out_eff = norm(rx, ry, rz);

                r3_a_out_eff = 1/(a_out_eff*a_out_eff*a_out_eff);
            }
        }
        /*---------------------------------------------------------------------------*\
            combinations
        \*---------------------------------------------------------------------------*/
        constexpr int Lin_idx = 0;
        constexpr int Lout_idx = 1;
        constexpr int S1_idx = 0;
        constexpr int S2_idx = 1;
        constexpr int S3_idx = 2;

        if constexpr (LL) {
            SL_coupling<false, Container, -4, Lout_idx>(args.args.Lin_Lout_coef * r3_a_out_eff, var, ddt);//evolve L1
            SL_coupling<false, Container, -3, Lout_idx>(args.args.Lin_Lout_coef * r3_a_out_eff, var, ddt);//evolve e1
        }

        if constexpr(spin_num(var) == 1) {
            SL_coupling<false, Container, S1_idx, Lin_idx>(args.S1_Lin_coef * r3_a_in_eff, var, ddt);
            if constexpr(SL) {
                SL_coupling<false, Container, S1_idx, Lout_idx>(args.S1_Lout_coef * r3_a_out_eff, var, ddt);
            }
        }

        if constexpr(spin_num(var) == 2) {
            SL_coupling<false, Container, S2_idx, Lin_idx>(args.S2_Lin_coef * r3_a_in_eff, var, ddt);
            if constexpr(SL) {
                SL_coupling<false, Container, S2_idx, Lout_idx>(args.S2_Lout_coef * r3_a_out_eff, var, ddt);
            }
        }

        if constexpr(spin_num(var) == 3) {
            SL_coupling<true, Container, S3_idx, L_out_idx>(args.S3_Lout_coef * r3_a_out_eff, var, ddt);
        }
    }

} // namespace secular
#endif
