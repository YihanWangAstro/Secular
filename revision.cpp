
#include <iostream>
#include "secular.h"
#include "SpaceHub/src/multi-thread/multi-thread.hpp"
#include "SpaceHub/src/tools/config-reader.hpp"
#include "SpaceHub/src/tools/timer.hpp"
#include "boost/numeric/odeint.hpp"
#include <iomanip>
#include <cstdlib>
#include <algorithm>

using namespace space::multiThread;

double INT_ERROR = 1e-13;

void get_line(std::fstream& is, std::string& str) {
    std::getline(is, str);
    if(!is)
      throw secular::StopFlag::eof;
}

void unpack_args_from_str(std::string  const& str, std::vector<double>& vec, bool DA, size_t spin_num) {
    std::stringstream is{str};

    size_t token_num = 22 + static_cast<size_t>(!DA) + spin_num*3;

    double tmp;

    vec.reserve(token_num);
    for(size_t i = 0 ;  i < token_num; ++i){
        is >> tmp;
        vec.emplace_back(tmp);
    }
}

auto resolve_sim_type(std::string const& line) {
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

    if(token_num>0) {
        id = std::stoi(line);
    } else{
        return std::make_tuple(id, double_average, 0u);
    }

    switch(token_num) {
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
            std::cout << "wrong input format!\n";
            throw secular::StopFlag::input_err;
    }
}


template<size_t spin_num, typename Observer>
void call_ode_int(bool DA, secular::Controler const&ctrl, secular::OrbitArgs const& init_args, double t_start, double t_end, Observer obsv){
    using Container = std::array<double, 12+spin_num*3>;
    using stepper_type = boost::numeric::odeint::runge_kutta_fehlberg78<Container>;

    Container init_cond;

    initilize_orbit_args(DA, spin_num, init_cond, init_args);

    secular::SecularConst const_parameters{ctrl, init_args.m1, init_args.m2, init_args.m3};

    auto func = secular::Dynamic_dispatch<Container>(ctrl, const_parameters);

    //auto func = secular::Static_dispatch<decltype(task.ctrl), decltype(args), Container>(task.ctrl, args);
    double ini_dt = 0.1 * secular::consts::year;

    boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::make_controlled(INT_ERROR, INT_ERROR, stepper_type()), func, init_cond, t_start, t_end, ini_dt, obsv);
}

constexpr size_t TASK_ID_OFFSET = 0;
constexpr size_t TRAJECTROY_OFFSET = 1;
constexpr size_t END_STAT_OFFSET = 2;
constexpr size_t END_TIME_OFFSET = 3;
constexpr size_t DT_OFFSET = 4;
constexpr size_t CTRL_OFFSET = 5;
constexpr size_t ARGS_OFFSET = 10;

struct Traj_args{
    Traj_args(std::string const& work_dir, size_t task_id, double _dt)
      : dt{_dt},
        t_output{0},
        file{work_dir +  "output_" + std::to_string(task_id) + ".txt", std::fstream::out} {
            file << std::setprecision(15);
        }

    void move_to_next_output(){ t_output += dt;}

    double dt;
    double t_output;
    std::fstream file;
};
/*
auto create_obser(std::shared_ptr<Traj_args>& traj_arg_ptr, std::string const& work_dir, std::vector<double> &input_args ) {
    bool is_traj = input_args[TRAJECTROY_OFFSET];
    size_t task_id = input_args[TASK_ID_OFFSET];
    double dt = input_args[DT_OFFSET];

    traj_arg_ptr = std::make_shared<Traj_args>(work_dir, task_id, dt);

    return [=](auto const& data, double t) {
        if(is_traj && t > traj_arg_ptr->t_output) {
            space::display(traj_arg_ptr->file, t, data, "\r\n");
            traj_arg_ptr->move_to_next_output();
        }
    };
}*/

void single_thread_job(std::string work_dir, ConcurrentFile input, size_t start_id, size_t end_id, ConcurrentFile output, ConcurrentFile log) {
    std::string entry;
    for (;;) {
        try{
            input.execute(get_line, entry);

            auto [task_id, DA, spin_num] = resolve_sim_type(entry);

            if(start_id<=task_id && task_id <= end_id){
                std::vector<double> v;

                unpack_args_from_str(entry, v, DA, spin_num);

                secular::Controler ctrl{v.begin() + CTRL_OFFSET, DA};

                log << get_log_title(task_id, DA, ctrl, spin_num);

                secular::OrbitArgs init_args{v.begin() + ARGS_OFFSET, DA, spin_num};

                double t_end = v[END_TIME_OFFSET];

                //space::display(std::cout, task_id, DA, spin_num, ctrl.Oct, ctrl.GR, ctrl.GW, ctrl.SL, ctrl.LL,"\n");

                bool is_traj = v[TRAJECTROY_OFFSET];

                size_t task_id = v[TASK_ID_OFFSET];

                double dt = v[DT_OFFSET];

                Traj_args traj_arg{work_dir, task_id, dt};

                auto observer = [&](auto const& data, double t) {
                    if(is_traj && t > traj_arg.t_output) {
                        traj_arg.file << t << ' ';
                        for(auto d : data) {
                            traj_arg.file << d << ' ';
                        }
                        traj_arg.file <<  "\r\n";
                        traj_arg.move_to_next_output();
                    }
                };

                if(spin_num == 0) {
                    call_ode_int<0>(DA, ctrl, init_args, 0.0, t_end, observer);
                } else if(spin_num == 1){
                    call_ode_int<1>(DA, ctrl, init_args, 0.0, t_end, observer);
                } else if(spin_num == 2){
                    call_ode_int<2>(DA, ctrl, init_args, 0.0, t_end, observer);
                } else if(spin_num == 3){
                    call_ode_int<3>(DA, ctrl, init_args, 0.0, t_end, observer);
                }
            }
        } catch (secular::StopFlag flag) {
            if(flag == secular::StopFlag::eof)
                break;
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
        exit(1);
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
