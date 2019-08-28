
#include <iostream>
#include "secular.h"
#include "tools.h"
#include "SpaceHub/src/multi-thread/multi-thread.hpp"
#include "SpaceHub/src/tools/config-reader.hpp"
#include "SpaceHub/src/tools/timer.hpp"
#include "boost/numeric/odeint.hpp"
#include <iomanip>
#include <cstdlib>
#include <algorithm>

using namespace space::multiThread;

double INT_ERROR = 1e-13;

bool get_line(std::fstream &is, std::string &str) {
    std::getline(is, str);
    if (!is)
        return false;
    else
        return true;
}

void unpack_args_from_str(std::string const &str, std::vector<double> &vec, bool DA, size_t spin_num) {
    std::stringstream is{str};

    size_t token_num = 22 + static_cast<size_t>(!DA) + spin_num * 3;

    double tmp;

    vec.reserve(token_num);
    for (size_t i = 0; i < token_num; ++i) {
        is >> tmp;
        vec.emplace_back(tmp);
    }
}

auto resolve_sim_type(std::string const &line) {
    bool in_space = true;
    size_t token_num = 0;
    for (auto c : line) {
        if (std::isspace(c)) {
            in_space = true;
        } else if (in_space) {
            token_num++;
            in_space = false;
        }
    }
    constexpr static bool single_average{false};
    constexpr static bool double_average{true};

    size_t id = 0;

    if (token_num > 0) {
        id = std::stoi(line);
    } else {
        return std::make_tuple(id, double_average, 0u);
    }

    switch (token_num) {
        case 22 :
            return std::make_tuple(id, double_average, 0u);
        case 23 :
            return std::make_tuple(id, single_average, 0u);
        case 25 :
            return std::make_tuple(id, double_average, 1u);
        case 26 :
            return std::make_tuple(id, single_average, 1u);
        case 28 :
            return std::make_tuple(id, double_average, 2u);
        case 29 :
            return std::make_tuple(id, single_average, 2u);
        case 31 :
            return std::make_tuple(id, double_average, 3u);
        case 32 :
            return std::make_tuple(id, single_average, 3u);
        default :
            return std::make_tuple(0lu, single_average, 0u);
    }
}

bool is_on(double x) {
    return x > 5e-15;
}
struct Stream_observer
{
    Stream_observer(std::ostream &out, double dt) : dt_{dt}, t_out_{0.0}, f_out_{out}, switch_{is_on(dt)} { }

    template<typename State>
    void operator()(State const&x , double t)
    {
        if(switch_ && t >= t_out_) {
              f_out_ << t << ' ' << x << "\r\n";
              t_out_ += dt_;
        }
    }
  private:
    double const dt_;
    double t_out_;
    std::ostream& f_out_;
    const bool switch_;
};


struct GW_Determinator
{
    GW_Determinator(double m1, double m2, double a_min) : a_min_{a_min}, a_coef_{0}  { }

    template<typename State>
    void operator()(State const&x , double t)
    {

    }
  private:
    double const a_min_;
    double const a_coef_;
    const bool switch_{true};
};

constexpr size_t TASK_ID_OFFSET = 0;
constexpr size_t TRAJECTROY_OFFSET = 1;
constexpr size_t GW_10HZ_OFFSET = 2;
constexpr size_t END_TIME_OFFSET = 3;
constexpr size_t DT_OFFSET = 4;
constexpr size_t CTRL_OFFSET = 5;
constexpr size_t ARGS_OFFSET = 10;

enum class ReturnFlag {
    shrink, eof, input_err, max_iter, finish
};

template<size_t spin_num>
auto call_ode_int(std::string work_dir, bool DA, secular::Controler const &ctrl, std::vector<double> const &init_args) {
    using namespace boost::numeric::odeint;

    using Container = secular::SecularArray<spin_num>;
    //using stepper_type = boost::numeric::odeint::runge_kutta_fehlberg78<Container>;
    using stepper_type = bulirsch_stoer<Container>;

    auto  [task_id, is_traj, is_10hz, t_end, out_dt] = secular::cast_unpack<decltype(init_args.begin()), size_t, bool, bool, double, double>(init_args.begin());

    auto const [m1, m2, m3, a_in_init] = secular::unpack_args<4>(init_args.begin() + ARGS_OFFSET);

    std::fstream f_out;

    if (is_on(out_dt))
        f_out.open(work_dir + "output_" + std::to_string(task_id) + ".txt", std::fstream::out);

    Container data;

    initilize_orbit_args(DA, spin_num, data, init_args.begin() + ARGS_OFFSET);

    secular::SecularConst<spin_num> const_parameters{m1, m2, m3};

    double dt = 0.1 * secular::consts::year;

    double time = 0;

    stepper_type stepper{INT_ERROR, INT_ERROR};

    auto ode = secular::Dynamic_dispatch<Container>(ctrl, const_parameters);

    Stream_observer writer{f_out, out_dt};

    for( ;time <= t_end; ) {

        writer(data, t_start);

        constexpr size_t max_attempts = 500;
        controlled_step_result res = success;
        size_t trials = 0;
        do{
            res = stepper.try_step(ode, data, time, dt );
            trails++;
        } while((res == fail) && trails < max_attempts);
        if(trails == max_attempts){
            return ReturnFlag::max_iter;
        }
    }

    //integrate_adaptive(stepper_type{INT_ERROR, INT_ERROR}, func, init_cond, t_start, t_end, ini_dt, obsv);
    //STATIC_DISPATH(ctrl, const_parameters, integrate_adaptive(stepper_type{INT_ERROR, INT_ERROR}, func, init_cond, t_start, t_end, ini_dt, Stream_observer(f_out, out_dt));)

    //STATIC_DISPATH(ctrl, const_parameters, integrate_adaptive(make_controlled(INT_ERROR, INT_ERROR, stepper_type()), func, init_cond, t_start, t_end, ini_dt, obsv);)

    //integrate_adaptive(make_controlled(INT_ERROR, INT_ERROR, stepper_type()), func, init_cond, t_start, t_end, ini_dt, obsv);

    return ReturnFlag::finish;
}

/*
template<size_t spin_num>
void call_ode_int(std::string work_dir, bool DA, secular::Controler const &ctrl, std::vector<double> const &init_args) {
    using namespace boost::numeric::odeint;

    using Container = secular::SecularArray<spin_num>;
    //using stepper_type = boost::numeric::odeint::runge_kutta_fehlberg78<Container>;
    using stepper_type = bulirsch_stoer<Container>;

    auto  [task_id, is_traj, is_10hz, t_end, out_dt] = secular::cast_unpack<decltype(init_args.begin()), size_t, bool, bool, double, double>(init_args.begin());

    auto const [m1, m2, m3, a_in_init] = secular::unpack_args<4>(init_args.begin() + ARGS_OFFSET);

    std::fstream f_out;

    if (is_on(out_dt))
        f_out.open(work_dir + "output_" + std::to_string(task_id) + ".txt", std::fstream::out);

    Container init_cond;

    initilize_orbit_args(DA, spin_num, init_cond, init_args.begin() + ARGS_OFFSET);

    secular::SecularConst<spin_num> const_parameters{m1, m2, m3};

    double ini_dt = 0.1 * secular::consts::year;

    double t_start = 0;

    auto func = secular::Dynamic_dispatch<Container>(ctrl, const_parameters);

    integrate_adaptive(stepper_type{INT_ERROR, INT_ERROR}, func, init_cond, t_start, t_end, ini_dt, Stream_observer(f_out, out_dt));
    //STATIC_DISPATH(ctrl, const_parameters, integrate_adaptive(stepper_type{INT_ERROR, INT_ERROR}, func, init_cond, t_start, t_end, ini_dt, Stream_observer(f_out, out_dt));)

    //STATIC_DISPATH(ctrl, const_parameters, integrate_adaptive(make_controlled(INT_ERROR, INT_ERROR, stepper_type()), func, init_cond, t_start, t_end, ini_dt, obsv);)

    //integrate_adaptive(make_controlled(INT_ERROR, INT_ERROR, stepper_type()), func, init_cond, t_start, t_end, ini_dt, obsv);
}*/

void single_thread_job(std::string work_dir, ConcurrentFile input, size_t start_id, size_t end_id, ConcurrentFile output, ConcurrentFile log) {
    std::string entry;
    for (;input.execute(get_line, entry);) {
            auto[task_id, DA, spin_num] = resolve_sim_type(entry);

            if (start_id <= task_id && task_id <= end_id) {
                std::vector<double> v;

                unpack_args_from_str(entry, v, DA, spin_num);

                secular::Controler ctrl{v.begin() + CTRL_OFFSET, DA};

                log << secular::get_log_title(task_id, DA, ctrl, spin_num) + "\r\n";

                if (spin_num == 0) {
                    call_ode_int<0>(work_dir, DA, ctrl, v);
                } else if (spin_num == 1) {
                    call_ode_int<1>(work_dir, DA, ctrl, v);
                } else if (spin_num == 2) {
                    call_ode_int<2>(work_dir, DA, ctrl, v);
                } else if (spin_num == 3) {
                    call_ode_int<3>(work_dir, DA, ctrl, v);
                }
            }
    }
}

int main(int argc, char **argv) {
    std::string input_file_name;
    std::string work_dir;
    size_t start_task_id, end_task_id;

    space::tools::read_command_line(argc, argv, input_file_name, start_task_id, end_task_id, work_dir);

    const int dir_err = system(("mkdir -p " + work_dir).c_str());
    if (dir_err == -1) {
        std::cout << "Error creating directory!\n";
        return 0;
    }

    work_dir += "/";

    auto input_file = make_thread_safe_fstream(input_file_name, std::fstream::in);

    auto output_file = make_thread_safe_fstream(work_dir + "statistics.txt", std::fstream::out);

    auto log_file = make_thread_safe_fstream(work_dir + "log.txt", std::fstream::out);

    std::cout << std::setprecision(15);

    space::tools::Timer timer;
    timer.start();
    space::multiThread::auto_multi_thread(single_thread_job, work_dir, input_file, start_task_id, end_task_id, output_file, log_file);
    std::cout << "\r\n Time:" << timer.get_time() << '\n';
    return 0;
}
