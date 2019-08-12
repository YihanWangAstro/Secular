#ifndef DESITTER_H
#define DESITTER_H

#include "tools.h"

namespace secular
{

inline double dsdt_coupling(double C, double a, double j_sqr){
    double r_a = 1.0/a;
    double r_a2 = r_a * r_a;
    double r_a4 = r_a2 * r_a2;
    double r_a5 = r_a4 * r_a;
    return C*sqrt(r_a5)/j_sqr;
}

inline auto self_deSitter(double Omega, double nx, double ny, double nz, double sx, double sy, double sz){
    auto const [cns_x, cns_y, cns_z] = cross(nx, ny, nz, sx, sy, sz);
    return std::make_tuple(Omega*cns_x, Omega*cns_y, Omega*cns_z);
}

template<typename Container>
constexpr size_t spin_num(Container const& d){
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

template<bool LL, bool SL, typename Args, typename Container>
inline void deSitter_precession(Args const& args, Container const& var, Container& ddt, double t){
        /*---------------------------------------------------------------------------*\
            mapping alias
        \*---------------------------------------------------------------------------*/
        const auto [L1x, L1y, L1z] = std::tie(var[0], var[1], var[2]);

        const auto [e1x, e1y, e1z] = std::tie(var[3], var[4], var[5]);
        /*---------------------------------------------------------------------------*\
            orbital parameters calculation
        \*---------------------------------------------------------------------------*/
        auto [e1_sqr, j1_sqr, j1, L1_norm, L_in, a_in] = calc_orbit_args(args.a_in_coef, L1x, L1y, L1z, e1x, e1y, e1z);
        /*---------------------------------------------------------------------------*\
            unit vectors
        \*---------------------------------------------------------------------------*/
        double n1x = L1x/L1_norm, n1y = L1y/L1_norm, n1z = L1z/L1_norm;  
        /*---------------------------------------------------------------------------*\
            combinations
        \*---------------------------------------------------------------------------*/

        if constexpr(spin_num(var) == 1){
            
            const auto [s1x, s1y, s1z] = std::tie(var[12], var[13], var[14]);

            auto [ds1x, ds1y, ds1z] = std::tie(ddt[12], ddt[13], ddt[14]);

            std::tie(ds1x, ds1y, ds1z) = self_deSitter(dsdt_coupling(args.S_1_L_in_coef, a_in, j1_sqr), n1x, n1y, n1z, s1x, s1y, s1z);
        } 
        
        if constexpr(spin_num(var) == 2) {
            const auto [s2x, s2y, s2z] = std::tie(var[15], var[16], var[17]);

            auto [ds2x, ds2y, ds2z] = std::tie(ddt[15], ddt[16], ddt[17]);

            std::tie(ds2x, ds2y, ds2z) = self_deSitter(dsdt_coupling(args.S_2_L_in_coef, a_in, j1_sqr), n1x, n1y, n1z, s2x, s2y, s2z);
        } 
        
        /*if constexpr(spin_num(var) == 3) {
            const auto [s3x, s3y, s3z] = std::tie(var[18], var[19], var[20]);

            auto [ds3x, ds3y, ds3z] = std::tie(ddt[18], ddt[19], ddt[20]);

            std::tie(ds3x, ds3y, ds3z) = self_deSitter(dsdt_coupling(args.S_3_L_out_coef, a_out, j2_sqr), n2x, n2y, n2z, s3x, s3y, s3z);
        }*/

        if constexpr(SL || LL) {
            /*---------------------------------------------------------------------------*\
                mapping alias
            \*---------------------------------------------------------------------------*/
            const auto [L2x, L2y, L2z] = std::tie(var[6], var[7], var[8]);

            const auto [e2x, e2y, e2z] = std::tie(var[9], var[10], var[11]);

            const auto [s1x, s1y, s1z] = std::tie(var[12], var[13], var[14]);

            auto [ds1x, ds1y, ds1z] = std::tie(ddt[12], ddt[13], ddt[14]);
            /*---------------------------------------------------------------------------*\
                orbital parameters calculation
            \*---------------------------------------------------------------------------*/
            auto [e2_sqr, j2_sqr, j2, L2_norm, L_out, a_out] = calc_orbit_args(args.a_out_coef, L2x, L2y, L2z, e2x, e2y, e2z);
            /*---------------------------------------------------------------------------*\
                unit vectors
            \*---------------------------------------------------------------------------*/
            double n2x = L2x/L2_norm, n2y = L2y/L2_norm, n2z = L2z/L2_norm;
            /*---------------------------------------------------------------------------*\
                combinations
            \*---------------------------------------------------------------------------*/
            /*if constexpr(var.size() >= 15){
                auto [ds1x_, ds1y_, ds1z_] = self_deSitter(dsdt_coupling(args.S_1_L_out_coef, a_out, j2_sqr), n2x, n2y, n2z, s1x, s1y, s1z);

                ds1x += ds1x_, ds1y += ds1y_, ds1z += ds1z_; 
            }

            if constexpr(var.size() >= 18){
                auto [ds2x_, ds2y_, ds2z_] = self_deSitter(dsdt_coupling(args.S_2_L_out_coef, a_out, j2_sqr), n2x, n2y, n2z, s2x, s2y, s2z);

                ds2x += ds2x_, ds2y += ds2y_, ds2z += ds2z_;
            }*/
        } 
}

} // namespace secular
#endif