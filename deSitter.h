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


    auto deSitter_increase(double Omega, double v1x, double v1y, double v1z, double v2x, double v2y, double v2z) {
        auto const[c_x, c_y, c_z] = cross(v1x, v1y, v1z, v2x, v2y, v2z);
        return std::make_tuple(Omega * c_x, Omega * c_y, Omega * c_z);
    }

    auto deSitter_e(double Omega, double v1x, double v1y, double v1z, double v2x, double v2y, double v2z)

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

        if constexpr (LL) {
          auto[dL1x, dL1y, dL1z] = std::tie(ddt[0], ddt[1], ddt[2]);

          auto[de1x, de1y, de1z] = std::tie(ddt[3], ddt[4], ddt[5]);

          double Omega_Lin_Lout = args.args.Lin_Lout_coef * r3_a_out_eff;

          auto[dL1x_, dL1y_, dL1z_] = deSitter_increase(Omega_Lin_Lout, L2x, L2y, L2z, L1x, L1y, L1z);

          dL1x += dL1x_, dL1y += dL1y_, dL1z += dL1z_;

          auto[de1x_, de1y_, de1z_] = deSitter_increase(Omega_Lin_Lout, L2x, L2y, L2z, e1x, e1y, e1z);

          de1x += de1x_, de1y += de1y_, de1z += de1z_;
        }

        /*---------------------------------------------------------------------------*\
            combinations
        \*---------------------------------------------------------------------------*/
        if constexpr(spin_num(var) == 1) {
            const auto[s1x, s1y, s1z] = std::tie(var[12], var[13], var[14]);

            auto[ds1x, ds1y, ds1z] = std::tie(ddt[12], ddt[13], ddt[14]);

            double Omega_S1_Lin = args.S1_Lin_coef * r3_a_in_eff;

            std::tie(ds1x, ds1y, ds1z) = deSitter_increase(Omega_S1_Lin, L1x, L1y, L1z, s1x, s1y, s1z);

            if constexpr(SL) {
                double Omega_S1_Lout = args.S1_Lout_coef * r3_a_out_eff;

                auto[ds1x_, ds1y_, ds1z_] = deSitter_increase(Omega_S1_Lout, L2x, L2y, L2z, s1x, s1y, s1z);

                ds1x += ds1x_, ds1y += ds1y_, ds1z += ds1z_;
            }
        }

        if constexpr(spin_num(var) == 2) {
            const auto[s2x, s2y, s2z] = std::tie(var[15], var[16], var[17]);

            auto[ds2x, ds2y, ds2z] = std::tie(ddt[15], ddt[16], ddt[17]);

            double Omega_S2_Lin = args.S2_Lin_coef * r3_a_in_eff;

            std::tie(ds2x, ds2y, ds2z) = deSitter_increase(Omega_S2_Lin, L1x, L1y, L1z, s2x, s2y, s2z);

            if constexpr(SL) {
                double Omega_S2_Lout = args.S2_Lout_coef * r3_a_out_eff;

                auto[ds2x_, ds2y_, ds2z_] = deSitter_increase(Omega_S2_Lout, L2x, L2y, L2z, s2x, s2y, s2z);

                ds2x += ds2x_, ds2y += ds2y_, ds2z += ds2z_;
            }
        }

        if constexpr(spin_num(var) == 3) {
            const auto[s3x, s3y, s3z] = std::tie(var[18], var[19], var[20]);

            auto[ds3x, ds3y, ds3z] = std::tie(ddt[18], ddt[19], ddt[20]);

            auto[dL2x, dL2y, dL2z] = std::tie(ddt[6], ddt[7], ddt[8]);

            auto[de2x, de2y, de2z] = std::tie(ddt[9], ddt[10], ddt[11]);

            double Omega_S3_Lout = args.S3_Lout_coef * r3_a_out_eff;

            std::tie(ds3x, ds3y, ds3z) = deSitter_increase(Omega_S3_Lout, L2x, L2y, L2z, s3x, s3y, s3z);

            
        }
    }

} // namespace secular
#endif
