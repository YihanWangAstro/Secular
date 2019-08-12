#ifndef RELATIVISTIC_H
#define RELATIVISTIC_H

#include "tools.h"

namespace secular
{

template<typename Args, typename Container>
inline void GR_precession(Args const& args, Container const& var, Container& ddt, double t){
        /*---------------------------------------------------------------------------*\
            mapping alias
        \*---------------------------------------------------------------------------*/
        const auto [L1x, L1y, L1z] = std::tie(var[0], var[1], var[2]);

        const auto [e1x, e1y, e1z] = std::tie(var[3], var[4], var[5]);

        auto [de1x, de1y, de1z] = std::tie(ddt[3], ddt[4], ddt[5]);
        /*---------------------------------------------------------------------------*\
            orbital parameters calculation
        \*---------------------------------------------------------------------------*/
        auto [e1_sqr, j1_sqr, j1, L1_norm, L_in, a_in] = calc_orbit_args(args.a_in_coef, L1x, L1y, L1z, e1x, e1y, e1z);
        /*---------------------------------------------------------------------------*\
            unit vectors
        \*---------------------------------------------------------------------------*/
        double n1x = L1x/L1_norm, n1y = L1y/L1_norm, n1z = L1z/L1_norm;
        /*---------------------------------------------------------------------------*\
            cross production
        \*---------------------------------------------------------------------------*/
        auto const [cn1e1_x, cn1e1_y, cn1e1_z] = cross(n1x, n1y, n1z, e1x, e1y, e1z);
        /*---------------------------------------------------------------------------*\
            combinations
        \*---------------------------------------------------------------------------*/
        double r_a = 1.0 / a_in;

        double GR_coef = args.GR_coef * sqrt(r_a * r_a * r_a * r_a * r_a)/ j1_sqr;

        de1x += GR_coef * cn1e1_x;

        de1y += GR_coef * cn1e1_y;

        de1z += GR_coef * cn1e1_z;
}

template<typename Args, typename Container>
inline void GW_radiation(Args const& args, Container const& var, Container& ddt, double t){
        /*---------------------------------------------------------------------------*\
            mapping alias
        \*---------------------------------------------------------------------------*/
        const auto [L1x, L1y, L1z] = std::tie(var[0], var[1], var[2]);

        const auto [e1x, e1y, e1z] = std::tie(var[3], var[4], var[5]);

        auto [dL1x, dL1y, dL1z] = std::tie(ddt[0], ddt[1], ddt[2]);

        auto [de1x, de1y, de1z] = std::tie(ddt[3], ddt[4], ddt[5]);
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
        double r_a = 1.0 / a_in;

        double r_a2 = r_a * r_a;

        double r_a4 = r_a2 * r_a2;

        double r_a7 = r_a4 * r_a2 * r_a;

        double r_j2 = 1.0 / j1_sqr;

        double r_j4 = r_j2 * r_j2;

        double r_j5 = r_j4 / j1;

        double GW_L_coef = args.GW_L_coef * sqrt(r_a7) * r_j4 * (1 + 0.875 * e1_sqr);

        double GW_e_coef = args.GW_e_coef * r_a4 * r_j5  * (1 + 121.0/304 * e1_sqr);

        dL1x += GW_L_coef * n1x;

        dL1y += GW_L_coef * n1y;

        dL1z += GW_L_coef * n1z;

        de1x += GW_e_coef * e1x;

        de1y += GW_e_coef * e1y;

        de1z += GW_e_coef * e1z;
}

} // namespace secular
#endif