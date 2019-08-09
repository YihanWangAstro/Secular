#ifndef SECULAR_H
#define SECULAR_H

#include <iostream>
#include <array>
#include <tuple>
#include <fstream>
#include <memory>
#include <cmath>
#include <functional>
#include <iomanip>
#include <boost/numeric/odeint.hpp>
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
        double tmp{0};
        is >> tmp; t.id = static_cast<size_t>(tmp);
        is >> tmp; t.write_traj = static_cast<bool>(tmp);
        is >> tmp; t.write_end = static_cast<bool>(tmp);
        is >> t.end_time;
        is >> t.dt_out;
        is >> tmp; t.Oct = static_cast<bool>(tmp);
        is >> tmp; t.GR = static_cast<bool>(tmp);
        is >> tmp; t.GW = static_cast<bool>(tmp);
        is >> tmp; t.SL_out = static_cast<bool>(tmp);
        is >> tmp; t.LL_couple = static_cast<bool>(tmp);

        //sis >> t.id >> t.write_traj >> t.write_end >> t.end_time >> t.dt_out >> t.Oct >> t.GR >> t.GW >> t.SL_out >> t.LL_couple;
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

inline double norm2(double x, double y, double z){
    return x*x + y*y + z*z;
}

inline double norm(double x, double y, double z){
    return sqrt(norm2(x,y,z));
}

inline double dot(double x1, double y1, double z1, double x2, double y2, double z2) {
    return x1*x2 + y1*y2 + z1*z2;
}

inline void cross(double& x3, double& y3, double& z3, double x1, double y1, double z1, double x2, double y2, double z2){
    x3 =  
    y3 = 
    z3 =
}

template <typename Container, typename Args>
struct Equations
{
public:
    using Stepper = std::function<void(Args const &, Container const &, Container&, double t)>;
    using Observer = std::function<void(Args const &, Container const &, double t)>;

    Equations(Args const &_args) : args{_args} {}

    void operator()(Container const &x, Container &dxdt, double t)
    {
        dxdt.fill(0);

        args.update(x);

        for(auto const& f: func){
            f(args, x, dxdt, t);
        }
    }

    void operator()(Container const &x, double t)
    {
        for(auto const& g: obs){
            g(args, x, t);
        }
    }

private:
    Args args;
    std::vector<Stepper> func;
    std::vector<Observer> obs;
};



template<typename Args, typename Container>
void DA_quad_kozai(Args const& args, Container const& var, Container& ddt, double t){
    double& L1x{var[0]};
    double& L1y{var[1]};
    double& L1z{var[2]};

    double& e1x{var[3]};
    double& e1y{var[4]};
    double& e1z{var[5]};

    double& L2x{var[6]};
    double& L2y{var[7]};
    double& L2z{var[8]};

    double& e2x{var[9]};
    double& e2y{var[10]};
    double& e2z{var[11]};

    double L1_norm = norm(L1x, L1y, L1z);

    double L2_norm = norm(L2x, L2y, L2z);

    double de1e1 = norm2(e1x, e1y, e1z);

    double a_out_eff = args.a_out * args.j2;

    double quad_coef = 0.75 / t_k_quad(args.m1 + args.m2, args.m3, args.a_in, a_out_eff);
}

template<typename Args, typename Container, typename ...Func>
auto serilize(Func&&...func){
    return [=](Container const &x, Container &dxdt, double t){ 
        (func(),...);
    };
}


template <typename ... Args>
auto f(Args&& ... args){
    return [args = std::make_tuple(std::forward<Args>(args) ...)]()mutable{
        return std::apply([](auto&& ... args){
            // use args
        }, std::move(args));
    };
}

} // namespace secular
#endif