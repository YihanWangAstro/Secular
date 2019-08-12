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

    auto deSitter_increase(double coef, double a_eff, double v1x, double v1y, double v1z, double v2x, double v2y,
                           double v2z) {
        auto const[c_x, c_y, c_z] = cross(v1x, v1y, v1z, v2x, v2y, v2z);
        double Omega = coef / (a_eff * a_eff * a_eff);
        return std::make_tuple(Omega * c_x, Omega * c_y, Omega * c_z);
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

        double a_out_eff{0};

        if constexpr(LL || SL || spin_num(var) == 3) {//if out orbit is coupled, we need to calculate a_out_eff
            if constexpr(DA) {//double averaged
                L2x = var[6], L2y = var[7], L2z = var[8];

                const auto[e2x, e2y, e2z] = std::tie(var[9], var[10], var[11]);

                a_out_eff = calc_a_eff(args.a_out_coef, L2x, L2y, L2z, e2x, e2y, e2z);
            } else {//single averaged
                const auto[rx, ry, rz] = std::tie(var[6], var[7], var[8]);

                const auto[vx, vy, vz] = std::tie(var[9], var[10], var[11]);

                std::tie(L2x, L2y, L2z) = cross_with_coef(args.mu2, rx, ry, rz, vx, vy, vz);

                a_out_eff = norm(rx, ry, rz);
            }
        }

        /*---------------------------------------------------------------------------*\
            combinations
        \*---------------------------------------------------------------------------*/

        if constexpr(spin_num(var) == 1) {

            const auto[s1x, s1y, s1z] = std::tie(var[12], var[13], var[14]);

            auto[ds1x, ds1y, ds1z] = std::tie(ddt[12], ddt[13], ddt[14]);

            std::tie(ds1x, ds1y, ds1z) = deSitter_increase(args.S_1_L_in_coef, a_in_eff, L1x, L1y, L1z, s1x, s1y, s1z);

            if constexpr(SL) {
                auto[ds1x_, ds1y_, ds1z_] = deSitter_increase(args.S_1_L_out_coef, a_out_eff, L2x, L2y, L2z, s1x, s1y,
                                                              s1z);

                ds1x += ds1x_, ds1y += ds1y_, ds1z += ds1z_;
            }
        }

        if constexpr(spin_num(var) == 2) {
            const auto[s2x, s2y, s2z] = std::tie(var[15], var[16], var[17]);

            auto[ds2x, ds2y, ds2z] = std::tie(ddt[15], ddt[16], ddt[17]);

            std::tie(ds2x, ds2y, ds2z) = deSitter_increase(args.S_2_L_in_coef, a_in_eff, L1x, L1y, L1z, s2x, s2y, s2z);

            if constexpr(SL) {
                auto[ds2x_, ds2y_, ds2z_] = deSitter_increase(args.S_2_L_out_coef, a_out_eff, L2x, L2y, L2z, s2x, s2y,
                                                              s2z);

                ds2x += ds2x_, ds2y += ds2y_, ds2z += ds2z_;
            }
        }

        if constexpr(spin_num(var) == 3) {
            const auto[s3x, s3y, s3z] = std::tie(var[18], var[19], var[20]);

            auto[ds3x, ds3y, ds3z] = std::tie(ddt[18], ddt[19], ddt[20]);

            std::tie(ds3x, ds3y, ds3z) = deSitter_increase(args.S_3_L_out_coef, a_out_eff, L2x, L2y, L2z, s3x, s3y,
                                                           s3z);
        }
    }

} // namespace secular
#endif