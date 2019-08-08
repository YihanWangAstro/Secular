#ifndef SECULAR_H
#define SECULAR_H

#include <iostream>
#include <array>
#include <tuple>
#include <fstream>
#include <memory>
#include <cmath>
#include <iomanip>
#include <boost/numeric/odeint.hpp>
#include "vector3d.h"
#include "SpaceHub/src/multi-thread/multi-thread.hpp"
#include "SpaceHub/src/dev-tools.hpp"

namespace secular
{

constexpr double pi = 3.14159265358979323;
constexpr double G = 4 * pi * pi;
constexpr double C = 6.32397263e4;
constexpr double r_G_sqrt = 1.0 / (2 * pi);

constexpr double year = 1;

double calc_angular_mom(double m_in, double m_out, double a)
{
    double mu = m_in * m_out / (m_in + m_out);
    return mu * sqrt(G * (m_in + m_out) * a);
}

Vec3d unit_j(double i, double Omega)
{
    double sini = sin(i);
    double cosi = cos(i);
    double sinO = sin(Omega);
    double cosO = cos(Omega);
    return Vec3d{sini * sinO, -sini * cosO, cosi};
}

Vec3d unit_e(double i, double omega, double Omega)
{
    double sini = sin(i);
    double cosi = cos(i);
    double sinO = sin(Omega);
    double cosO = cos(Omega);
    double sino = sin(omega);
    double coso = cos(omega);

    return Vec3d{coso * cosO - sino * cosi * sinO, coso * sinO + sino * cosi * cosO, sino * sini};
}

double t_k_quad(double m_in, double m_out, double a_in, double a_out_eff)
{
    double ratio = a_out_eff/ sqrt(a_in);
    return r_G_sqrt * sqrt(m_in) / m_out * ratio * ratio * ratio;
}

double normed_oct_epsilon(double m1, double m2, double a_in, double a_out, double c_out_sqr) {
    return fabs(m1-m2)/(m1+m2)*a_in/a_out/c_out_sqr;
}

template <typename T, size_t len>
std::ostream &operator<<(std::ostream &os, std::array<T, len> const &arr)
{
    for (auto const &c : arr)
    {
        os << c << ' ';
    }
    return os;
}

template <size_t SpinNum_>
struct OrbitArgs
{
    static constexpr size_t SpinNum{SpinNum_};
    double m1;
    double m2;
    double m3;
    double a_in;
    double a_out;
    double e_in;
    double e_out;
    double omega_in;
    double omega_out;
    double Omega_in;
    double Omega_out;
    double i_in;
    double i_out;
    std::array<Vec3d, SpinNum> s;

    friend std::istream &operator>>(std::istream &is, OrbitArgs &t)
    {
        is >> t.m1 >> t.m2 >> t.m3 >> t.a_in >> t.a_out >> t.e_in >> t.e_out >> t.omega_in >> t.omega_out >> t.Omega_in >> t.i_in >> t.i_out;
        t.Omega_out = t.Omega_in - 180.0;
        for (auto &ss : t.s)
        {
            is >> ss;
        }
        return is;
    }
    friend std::ostream &operator<<(std::ostream &os, OrbitArgs const &t)
    {
        os << t.m1 << ' ' << t.m2 << ' ' << t.m3 << ' ' << t.a_in << ' ' << t.a_out << ' ' << t.e_in << ' ' << t.e_out << ' ' << t.omega_in << ' ' << t.omega_out << ' ' << t.Omega_in << ' ' << t.Omega_out << ' ' << t.i_in << ' ' << t.i_out;
        for (auto const &ss : t.s)
        {
            os << ' ' << ss;
        }
        return os;
    }
};

template <typename Obt>
void deg_to_rad(Obt &args)
{
    constexpr double rad = pi / 180.0;
    args.omega_in *= rad;
    args.omega_out *= rad;
    args.Omega_in *= rad;
    args.Omega_out *= rad;
    args.i_in *= rad;
    args.i_out *= rad;
}

struct Controler
{
    size_t id;
    bool write_traj;
    bool write_end;
    double end_time;
    double dt_out;
    bool Oct;
    bool GR;
    bool GW;
    bool SL_out;
    bool LL_couple;

    friend std::istream &operator>>(std::istream &is, Controler &t)
    {
        is >> t.id >> t.write_traj >> t.write_end >> t.end_time >> t.dt_out;
        double tmp{0};
        is >> tmp; t.Oct = static_cast<bool>(tmp);
        is >> tmp; t.GR = static_cast<bool>(tmp);
        is >> tmp; t.GW = static_cast<bool>(tmp);
        is >> tmp; t.SL_out = static_cast<bool>(tmp);
        is >> tmp; t.LL_couple = static_cast<bool>(tmp);
        return is;
    }

    friend std::ostream &operator<<(std::ostream &os, Controler const &t)
    {
        os << t.id << ' ' << t.write_traj << ' ' << t.write_end << ' ' << t.end_time << ' ' << t.dt_out << ' ' << t.Oct << ' ' << t.GR << ' ' << t.GW << ' ' << t.SL_out << ' ' << t.LL_couple;
        return os;
    }
};

template <size_t spin_num>
struct Task
{
    Controler ctrl;
    OrbitArgs<spin_num> obt_args;
    friend std::istream &operator>>(std::istream &is, Task &t)
    {
        is >> t.ctrl >> t.obt_args;
        return is;
    }
    friend std::ostream &operator<<(std::ostream &os, Task const &t)
    {
        os << t.ctrl << ' ' << t.obt_args;
        return os;
    }
};

enum class StopFlag
{
    shrink
};

template <size_t SpinNum>
struct State
{
    using Container = std::array<double, 12 + 3 * SpinNum>;

    State() = default;
    State(State const &) = default;
    State(State &&) = default;

    State(Container const &d)
        : L1{d[0], d[1], d[2]},
          e1{d[3], d[4], d[5]},
          L2{d[6], d[7], d[8]},
          e2{d[9], d[10], d[11]}
    {
        if (SpinNum > 0)
        {
            s[0].x = d[12], s[0].y = d[13], s[0].z = d[14];
        }
        if (SpinNum > 1)
        {
            s[1].x = d[15], s[1].y = d[16], s[1].z = d[17];
        }
        if (SpinNum > 2)
        {
            s[2].x = d[18], s[2].y = d[19], s[2].z = d[20];
        }
    }

    State &operator=(State const &d)
    {
        L1 = d.L1;
        e1 = d.e1;
        L2 = d.L2;
        e2 = d.e2;
        s = d.s;
        return *this;
    }

    State &operator=(Container const &d)
    {
        L1.x = d[0], L1.y = d[1], L1.z = d[2];
        e1.x = d[3], e1.y = d[4], e1.z = d[5];
        L2.x = d[6], L2.y = d[7], L2.z = d[8];
        e2.x = d[9], e2.y = d[10], e2.z = d[11];
        if (SpinNum > 0)
        {
            s[0].x = d[12], s[0].y = d[13], s[0].z = d[14];
        }
        if (SpinNum > 1)
        {
            s[1].x = d[15], s[1].y = d[16], s[1].z = d[17];
        }
        if (SpinNum > 2)
        {
            s[2].x = d[18], s[2].y = d[19], s[2].z = d[20];
        }
        return *this;
    }

    friend void to_container(Container &x, State const &d)
    {
        x[0] = d.L1.x, x[1] = d.L1.y, x[2] = d.L1.z;
        x[3] = d.e1.x, x[4] = d.e1.y, x[5] = d.e1.z;
        x[6] = d.L2.x, x[7] = d.L2.y, x[8] = d.L2.z;
        x[9] = d.e2.x, x[10] = d.e2.y, x[11] = d.e2.z;
        if (SpinNum > 0)
        {
            x[12] = d.s[0].x, x[13] = d.s[0].y, x[14] = d.s[0].z;
        }
        if (SpinNum > 1)
        {
            x[15] = d.s[1].x, x[16] = d.s[1].y, x[17] = d.s[1].z;
        }
        if (SpinNum > 2)
        {
            x[18] = d.s[2].x, x[19] = d.s[2].y, x[20] = d.s[2].z;
        }
    }

    Vec3d L1;
    Vec3d e1;
    union{
        struct {
            Vec3d L2;
            Vec3d e2;
        };
        struct{
            Vec3d r_out;
            Vec3d v_out;
        };
    };
    
    std::array<Vec3d, SpinNum> s;
};

template <size_t SpinNum_, typename Ctrl>
struct Secular
{
public:
    //Typemember
    constexpr static size_t SpinNum{SpinNum_};
    using States = State<SpinNum>;
    using Container = typename States::Container;

    Secular(std::string &work_dir, space::multiThread::ConcurrentFile output, OrbitArgs<SpinNum> const &orbit, Ctrl &ctr)
        : m1{orbit.m1},
          m2{orbit.m2},
          m3{orbit.m3},
          a_in_init{orbit.a_in},
          mu1{m1 * m2 / (m1 + m2)},
          mu2{(m1 + m2) * m3 / (m1 + m2 + m3)},
          ctrl{ctr},
          f_stat_{output}
    {
        using namespace secular;

        a_in_coef = 1 / (G * (m1 + m2)) / mu1 / mu1;
        a_out_coef = 1 / (G * (m1 + m2 + m3)) / mu2 / mu2;

        double m12 = m1 + m2;
        double C5 = C * C * C * C * C;
        gw_L_coef = -6.4 * pow(G, 3.5) * mu1 * mu1 * pow(m12, 2.5) / C5;
        gw_e_coef = -304.0 / 15 * G * G * G * mu1 * m12 * m12 / C5;

        gr_coef = 3 * pow(G * m12, 1.5) / C / C;

        auto deSitter = [](double m_self, double m_other) {
            double m_tot = m_self + m_other;
            double mu = m_self * m_other / (m_self + m_other);

            return 1.5 * pow(G, 1.5) * sqrt(m_tot) * (m_other + mu / 3) / (C * C);
        };

        S_1_L_in_coef = deSitter(m1, m2);
        S_2_L_in_coef = deSitter(m2, m1);

        S_1_L_out_coef = deSitter(m1 + m2, m3);
        S_2_L_out_coef = S_1_L_out_coef;

        S_3_L_out_coef = deSitter(m3, m1 + m2);

        if (ctr.write_traj)
        {
            fout_ = std::make_shared<std::fstream>(work_dir + "trajectory_" + std::to_string(ctr.id) + ".txt", std::fstream::out);

            if (!fout_->is_open())
            {
                std::cout << "Fail to open the file!\n";
                exit(0);
            }
            else
            {
                (*fout_) << std::fixed << std::setprecision(14);
            }
        }

        States args;

        args.L1 = secular::calc_angular_mom(orbit.m1, orbit.m2, orbit.a_in) * sqrt(1 - orbit.e_in * orbit.e_in) * secular::unit_j(orbit.i_in, orbit.Omega_in);
        args.L2 = secular::calc_angular_mom(orbit.m1 + orbit.m2, orbit.m3, orbit.a_out) * sqrt(1 - orbit.e_out * orbit.e_out) * secular::unit_j(orbit.i_out, orbit.Omega_out);
        args.e1 = orbit.e_in * secular::unit_e(orbit.i_in, orbit.omega_in, orbit.Omega_in);
        args.e2 = orbit.e_out * secular::unit_e(orbit.i_out, orbit.omega_out, orbit.Omega_out);

        args.s = orbit.s;

        to_container(this->initial_conds, args);
    }

    void operator()(Container const &x, Container &dxdt, double t)
    {
        States args{x};
        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        double L1_norm = norm(args.L1);

        double L2_norm = norm(args.L2);

        Vec3d n1 = args.L1 / L1_norm;

        Vec3d n2 = args.L2 / L2_norm;

        double de1e1 = norm2(args.e1);

        double de2e2 = norm2(args.e2);

        double dn1n2 = dot(n1, n2);

        double de1n2 = dot(args.e1, n2);

        Vec3d cn1n2 = cross(n1, n2);

        Vec3d cn1e1 = cross(n1, args.e1);

        Vec3d ce1n2 = cross(args.e1, n2);

        Vec3d ce2n1 = cross(args.e2, n1);

        Vec3d ce1e2 = cross(args.e1, args.e2);

        Vec3d cn2e2 = cross(n2, args.e2);

        double c_in_sqr = 1 - de1e1;

        double c_out_sqr = 1 - de2e2;

        double c_in = sqrt(c_in_sqr);

        double c_out = sqrt(c_out_sqr);

        double L_in = L1_norm / c_in;

        double L_out = L2_norm / c_out;

        double a_in = calc_a_in(L_in);

        double a_out = calc_a_out(L_out);
        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        States d_args;

        //quadrupole LK
        double t_k = t_k_quad(m1 + m2, m3, a_in, a_out*c_out);

        double coef = 0.75 / t_k;
        
        d_args.L1 = coef * L_in * ((c_in_sqr * dn1n2) * cn1n2 - (5 * de1n2) * ce1n2);

        d_args.L2 = -d_args.L1;

        d_args.e1 = coef * c_in * (dn1n2 * ce1n2 + 2 * cn1e1 - 5 * de1n2 * cn1n2);

        d_args.e2 = coef * (L_in / L_out) / c_out * ( (5 * de1n2) * ce1e2 + (c_in_sqr * dn1n2) * ce2n1 - (0.5 - 3 * de1e1 + 12.5 * de1n2 * de1n2 - 2.5 * c_in_sqr * dn1n2 * dn1n2) * cn2e2);

        //Octupole 

        if(ctrl.Oct) {
            double oct_coef = -75.0/64 * normed_oct_epsilon(m1, m2, a_in, a_out, c_out_sqr)/t_k;

            double de1e2 = dot(args.e1, args.e2);

            double dn1e2 = dot(n1, d_args.e2);

            double C1 = (1.6 * de1e1 - 0.2 - 7 * de1n2 * de1n2 + c_in_sqr * dn1n2 * dn1n2);

            double C2 = de1e2 * dn1n2 + de1n2 * dn1e2;

            double C3 = c_in_sqr * dn1n2 * dn1e2 - 7 * de1n2 * de1e2;

            double C4 = de1n2 * dn1n2;
                  
            double _2c_in_sqr = 2 * c_in_sqr;

            Vec3d  oct_dL1dt = oct_coef * L_in * (_2c_in_sqr * ( C2 * cn1n2  - C4 * ce2n1 ) + 2 * C3 * ce1n2  + C1 * ce1e2);

            d_args.L1 += oct_dL1dt;

            //d_args.L2 -= oct_dL1dt;

            //d_args.e1 += oct_coef * 2 * c_in * (C4 * ce1e2 - (0.5 * C1) * ce2n1 + C2 * ce1n2 + C3 * cn1n2 + (1.6 * de1e2) * cn1e1);

            //d_args.e2 += oct_coef * (L_in / L_out) / c_out * (_2c_in_sqr * ( C2 * ce2n1 - (c_out_sqr * C4) * cn1n2) - (c_out_sqr * C1) * ce1n2 - 2 * C3 * ce1e2 + ((0.4 - 3.2*de1e1)*de1e2 + 14*dn1e2* C4 *c_in_sqr + 7*de1e2*C1) * cn2e2 );*/
        }

        //GW radiation
        if (ctrl.GW)
        {
            d_args.e1 += calc_gw_dedt(args.e1, a_in, c_in);

            d_args.L1 += calc_gw_dLdt(n1, a_in, c_in);
        }

        //GR
        if (ctrl.GR)
        {
            d_args.e1 += calc_gr_dedt(cn1e1, a_in, c_in);
        }

        //LL
        if (ctrl.LL_couple)
        {
            d_args.L1 += calc_coupling_dsdt(args.L1, n2, a_out, c_out, S_1_L_out_coef);

            d_args.e1 += calc_coupling_dsdt(args.e1, n2, a_out, c_out, S_1_L_out_coef); //note! this line calls the correct function. maybe confused
        }

        //Spin
        if (SpinNum > 0)
        {
            d_args.s[0] = calc_coupling_dsdt(args.s[0], n1, a_in, c_in, S_1_L_in_coef);
            if (ctrl.SL_out)
            {
                d_args.s[0] += calc_coupling_dsdt(args.s[0], n2, a_out, c_out, S_1_L_out_coef);
            }
        }
        if (SpinNum > 1)
        {
            d_args.s[1] = calc_coupling_dsdt(args.s[1], n1, a_in, c_in, S_2_L_in_coef);
            if (ctrl.SL_out)
            {
                d_args.s[1] += calc_coupling_dsdt(args.s[1], n2, a_out, c_out, S_2_L_out_coef);
            }
        }
        if (SpinNum > 2)
        {
            d_args.s[2] = calc_coupling_dsdt(args.s[2], n2, a_out, c_out, S_3_L_out_coef);

            d_args.L2 -= d_args.s[2];

            d_args.e2 += calc_coupling_dedt(args.s[2], n2, args.e2, L2_norm, a_out, c_out, S_3_L_out_coef);
        }

        to_container(dxdt, d_args);
    }

    void operator()(Container const &data, double t)
    {
        Vec3d L1{data[0], data[1], data[2]};
        Vec3d e1{data[3], data[4], data[5]};

        double e_sqr = norm2(e1);
        double L_sqr = norm2(L1);
        double a = a_in_coef * L_sqr / (1 - e_sqr);

        if (a < 0.001 * a_in_init)
        {
            if (ctrl.write_traj)
                (*fout_) << t << ' ' << data << std::endl;

            if (ctrl.write_end)
            {
                f_stat_ << PACK(ctrl.id, ' ', t, ' ', data, "\r\n");
                f_stat_.flush();
            }

            throw StopFlag::shrink;        }

        if (ctrl.write_traj && t >= out_time_)
        {
            (*fout_) << t << ' ';
            for (auto const &d : data)
            {
                (*fout_) << d << ' ';
            }
            (*fout_) << "\r\n";
            out_time_ += ctrl.dt_out;
        }
    }

private:
    inline double calc_a_in(double L)
    {
        return a_in_coef * L * L;
    }

    inline double calc_a_out(double L)
    {
        return a_out_coef * L * L;
    }

    inline Vec3d calc_gw_dLdt(Vec3d const &n1, double a, double c_in)
    {
        double e_sqr = 1 - c_in * c_in;
        double r_a = 1.0 / a;
        return (gw_L_coef * sqrt(r_a * r_a * r_a * r_a * r_a * r_a * r_a) / (c_in * c_in * c_in * c_in) * (1 + 0.875 * e_sqr)) * n1;
        //return (gw_L_coef * pow(a, -3.5) / (c_in * c_in * c_in * c_in) * (1 + 7.0/8 * e_sqr) )* n1;
    }

    inline Vec3d calc_gw_dedt(Vec3d const &e1, double a, double c_in)
    {
        double e_sqr = norm2(e1);
        double r_c = 1.0 / c_in;
        return (gw_e_coef * (1 + 121.0 / 304 * e_sqr) * (r_c * r_c * r_c * r_c * r_c) / (a * a * a * a)) * e1;
        //return (gw_e_coef * (1 + 121.0/304 *e_sqr) * pow(c_in, -5) / (a * a * a * a) ) * e1;
    }

    inline Vec3d calc_gr_dedt(Vec3d const &cn1e1, double a, double c_in)
    {
        double r_a = 1.0 / a;
        return (gr_coef * sqrt(r_a * r_a * r_a * r_a * r_a) / (c_in * c_in)) * cn1e1;
        //return (gr_coef * pow(a, -2.5)/ ( c_in * c_in ) ) * cn1e1;
    }

    inline Vec3d calc_coupling_dsdt(Vec3d const &S, Vec3d const &n, double a, double cir, double coef)
    {
        double r_a = 1.0 / a;
        return (coef * sqrt(r_a * r_a * r_a * r_a * r_a) / (cir * cir)) * cross(n, S);
        //return (coef * pow(a, -2.5) / (cir * cir)) * cross(n, S);
    }

    inline Vec3d calc_coupling_dedt(Vec3d const &S, Vec3d const &n, Vec3d const &e, double L, double a, double cir, double coef)
    {
        double r_a = 1.0 / a;
        return (coef * sqrt(r_a * r_a * r_a * r_a * r_a) / (cir * cir) / L) * (cross(S, e) - (3 * dot(n, S)) * cross(n, e));
        //return (coef * pow(a, -2.5) / (cir * cir) / L) * (cross(S, e) - (3 * dot(n, S)) * cross(n, e));
    }

public:
    Container initial_conds;
    const Ctrl ctrl;

private:
    const double m1;
    const double m2;
    const double m3;
    const double mu1;
    const double mu2;
    const double a_in_init;

    double a_in_coef;
    double a_out_coef;

    double gw_L_coef;
    double gw_e_coef;

    double gr_coef;

    double S_1_L_in_coef;
    double S_2_L_in_coef;

    double S_1_L_out_coef;
    double S_2_L_out_coef;
    double S_3_L_out_coef;

    std::shared_ptr<std::fstream> fout_{nullptr};
    space::multiThread::ConcurrentFile f_stat_;

    double out_time_{0};
};


template <size_t SpinNum_, typename Ctrl>
struct Secular_SA
{
public:
    //Typemember
    constexpr static size_t SpinNum{SpinNum_};
    using States = State<SpinNum>;
    using Container = typename States::Container;

    Secular_SA(std::string &work_dir, space::multiThread::ConcurrentFile output, OrbitArgs<SpinNum> const &orbit, Ctrl &ctr)
        : m1{orbit.m1},
          m2{orbit.m2},
          m3{orbit.m3},
          a_in_init{orbit.a_in},
          mu1{m1 * m2 / (m1 + m2)},
          mu2{(m1 + m2) * m3 / (m1 + m2 + m3)},
          ctrl{ctr},
          f_stat_{output}
    {
        using namespace secular;

        a_in_coef = 1 / (G * (m1 + m2)) / mu1 / mu1;

        double m12 = m1 + m2;
        double C5 = C * C * C * C * C;
        gw_L_coef = -6.4 * pow(G, 3.5) * mu1 * mu1 * pow(m12, 2.5) / C5;
        gw_e_coef = -304.0 / 15 * G * G * G * mu1 * m12 * m12 / C5;

        gr_coef = 3 * pow(G * m12, 1.5) / C / C;

        auto deSitter = [](double m_self, double m_other) {
            double m_tot = m_self + m_other;
            double mu = m_self * m_other / (m_self + m_other);

            return 1.5 * pow(G, 1.5) * sqrt(m_tot) * (m_other + mu / 3) / (C * C);
        };

        S_1_L_in_coef = deSitter(m1, m2);
        S_2_L_in_coef = deSitter(m2, m1);

        S_1_L_out_coef = deSitter(m1 + m2, m3);
        S_2_L_out_coef = S_1_L_out_coef;

        S_3_L_out_coef = deSitter(m3, m1 + m2);

        if (ctr.write_traj)
        {
            fout_ = std::make_shared<std::fstream>(work_dir + "trajectory_" + std::to_string(ctr.id) + ".txt", std::fstream::out);

            if (!fout_->is_open())
            {
                std::cout << "Fail to open the file!\n";
                exit(0);
            }
            else
            {
                (*fout_) << std::fixed << std::setprecision(14);
            }
        }

        States args;

        args.L1 = secular::calc_angular_mom(orbit.m1, orbit.m2, orbit.a_in) * sqrt(1 - orbit.e_in * orbit.e_in) * secular::unit_j(orbit.i_in, orbit.Omega_in);
        args.e1 = orbit.e_in * secular::unit_e(orbit.i_in, orbit.omega_in, orbit.Omega_in);


        //args.r_out = 0;
        //args.v_out = 0;

        args.s = orbit.s;

        to_container(this->initial_conds, args);
    }

    void operator()(Container const &x, Container &dxdt, double t)
    {
        States args{x};
        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        double L1_norm = norm(args.L1);

        double r_out_norm = norm(args.L2);

        Vec3d n1 = args.L1 / L1_norm;

        Vec3d r_out_hat = args.r_out / r_out_norm;

        double de1e1 = norm2(args.e1);

        double de1r = dot(args.e1, r_out_hat);

        double dn1r = dot(n1, r_out_hat);

        Vec3d  ce1r = cross(args.e1, r_out_hat);

        Vec3d  cn1r = cross(n1, r_out_hat);

        Vec3d  cn1e1 = cross(n1, args.e1);

        double c_in_sqr = 1 - de1e1;

        double c_in = sqrt(c_in_sqr);

        double L_in = L1_norm / c_in;

        double a_in = calc_a_in(L_in);

        
        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        States d_args;

        //quadrupole LK
        double t_k = t_k_quad(m1 + m2, m3, a_in, r_out_norm);

        double coef = 1.5 / t_k;
        
        d_args.L1 = coef * L_in * (5*de1r*ce1r - c_in_sqr*dn1r*cn1r);               

        d_args.e1 = coef * c_in * (5*de1r*cn1r - dn1r*ce1r - 2*cn1e1);

        //d_args.r_out = 0;

        //d_args.v_out = 0;

        //Octupole 

        if(ctrl.Oct) {
           

        //GW radiation
        if (ctrl.GW)
        {
            d_args.e1 += calc_gw_dedt(args.e1, a_in, c_in);

            d_args.L1 += calc_gw_dLdt(n1, a_in, c_in);
        }

        //GR
        if (ctrl.GR)
        {
            d_args.e1 += calc_gr_dedt(cn1e1, a_in, c_in);
        }

        //LL
        if (ctrl.LL_couple)
        {
            
        }

        //Spin
        if (SpinNum > 0)
        {
           
        }
        if (SpinNum > 1)
        {
            
        }
        if (SpinNum > 2)
        {
           
        }

        to_container(dxdt, d_args);
    }

    void operator()(Container const &data, double t)
    {
        Vec3d L1{data[0], data[1], data[2]};
        Vec3d e1{data[3], data[4], data[5]};

        double e_sqr = norm2(e1);
        double L_sqr = norm2(L1);
        double a = a_in_coef * L_sqr / (1 - e_sqr);

        if (a < 0.001 * a_in_init)
        {
            if (ctrl.write_traj)
                (*fout_) << t << ' ' << data << std::endl;

            if (ctrl.write_end)
            {
                f_stat_ << PACK(ctrl.id, ' ', t, ' ', data, "\r\n");
                f_stat_.flush();
            }

            throw StopFlag::shrink;        }

        if (ctrl.write_traj && t >= out_time_)
        {
            (*fout_) << t << ' ';
            for (auto const &d : data)
            {
                (*fout_) << d << ' ';
            }
            (*fout_) << "\r\n";
            out_time_ += ctrl.dt_out;
        }
    }

private:
    inline double calc_a_in(double L)
    {
        return a_in_coef * L * L;
    }

    inline Vec3d calc_gw_dLdt(Vec3d const &n1, double a, double c_in)
    {
        double e_sqr = 1 - c_in * c_in;
        double r_a = 1.0 / a;
        return (gw_L_coef * sqrt(r_a * r_a * r_a * r_a * r_a * r_a * r_a) / (c_in * c_in * c_in * c_in) * (1 + 0.875 * e_sqr)) * n1;
    }

    inline Vec3d calc_gw_dedt(Vec3d const &e1, double a, double c_in)
    {
        double e_sqr = norm2(e1);
        double r_c = 1.0 / c_in;
        return (gw_e_coef * (1 + 121.0 / 304 * e_sqr) * (r_c * r_c * r_c * r_c * r_c) / (a * a * a * a)) * e1;
        
    }

    inline Vec3d calc_gr_dedt(Vec3d const &cn1e1, double a, double c_in)
    {
        double r_a = 1.0 / a;
        return (gr_coef * sqrt(r_a * r_a * r_a * r_a * r_a) / (c_in * c_in)) * cn1e1;
    }

    inline Vec3d calc_coupling_dsdt(Vec3d const &S, Vec3d const &n, double a, double cir, double coef)
    {
        double r_a = 1.0 / a;
        return (coef * sqrt(r_a * r_a * r_a * r_a * r_a) / (cir * cir)) * cross(n, S);
    }

    inline Vec3d calc_coupling_dedt(Vec3d const &S, Vec3d const &n, Vec3d const &e, double L, double a, double cir, double coef)
    {
        double r_a = 1.0 / a;
        return (coef * sqrt(r_a * r_a * r_a * r_a * r_a) / (cir * cir) / L) * (cross(S, e) - (3 * dot(n, S)) * cross(n, e));
    }

public:
    Container initial_conds;
    const Ctrl ctrl;

private:
    const double m1;
    const double m2;
    const double m3;
    const double mu1;
    const double mu2;
    const double a_in_init;

    double a_in_coef;

    double gw_L_coef;
    double gw_e_coef;

    double gr_coef;

    double S_1_L_in_coef;
    double S_2_L_in_coef;

    double S_1_L_out_coef;
    double S_2_L_out_coef;
    double S_3_L_out_coef;

    std::shared_ptr<std::fstream> fout_{nullptr};
    space::multiThread::ConcurrentFile f_stat_;

    double out_time_{0};
};
} // namespace secular
#endif