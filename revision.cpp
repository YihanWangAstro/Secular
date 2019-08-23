
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

double int_error = 1e-13;
/*
size_t resolve_spin_num(std::string const &file_name, size_t base = 22) {
    std::fstream init(file_name, std::fstream::in);
    std::string line;
    std::getline(init, line);

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

    if (token_num == base) {
        return 0;
    } else if (token_num == base + 3) {
        return 1;
    } else if (token_num == base + 6) {
        return 2;
    } else if (token_num == base + 9) {
        return 3;
    } else {
        return 999;
    }
}


template<size_t spin_num>
void
single_thread_job(std::string work_dir, ConcurrentFile input, size_t start_id, size_t end_id, ConcurrentFile output,
                  ConcurrentFile log) {
    using Container = std::array<double, 12>;

    using stepper_type = boost::numeric::odeint::runge_kutta_fehlberg78<Container>;
    double ini_dt = 0.1 * secular::consts::year;
    secular::Task<spin_num> task;
    for (;;) {
        if (input >> task) {
            if (start_id <= task.ctrl.id && task.ctrl.id <= end_id) {
                secular::deg_to_rad(task.obt_args);

                std::array<double, 12> inits;

                if (task.ctrl.DA == true) {
                    initilize_DA(inits, task.obt_args);
                } else {
                    initilize_SA(inits, task.obt_args);
                }

                secular::SecularArg<decltype(task.ctrl), spin_num> args{task.ctrl, task.obt_args.m1, task.obt_args.m2, task.obt_args.m3};

                std::fstream fout_{work_dir + "trajectory_" + std::to_string(task.ctrl.id) + ".txt", std::fstream::out};

                fout_ << std::fixed << std::setprecision(14);

                double out_time_ = 0;


                auto observer = [&](auto &data, double t) {
                    if (t >= out_time_) {
                        (fout_) << t << ' ';
                        for (auto const &d : data) {
                            (fout_) << d << ' ';
                        }
                        (fout_) << "\r\n";
                        out_time_ += task.ctrl.dt_out;
                    }
                };
                //
                auto func = secular::Dynamic_dispatch<decltype(task.ctrl), decltype(args), Container>(task.ctrl, args);

                //auto func = secular::Static_dispatch<decltype(task.ctrl), decltype(args), Container>(task.ctrl, args);

                try {
                    boost::numeric::odeint::integrate_adaptive(
                            boost::numeric::odeint::make_controlled(int_error, int_error, stepper_type()), func, inits,
                            0.0, task.ctrl.end_time, ini_dt, observer);
                } catch (secular::StopFlag flag) {
                    if (flag == secular::StopFlag::shrink) {
                        //
                    }
                }
            }
        } else {
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

    size_t spin_num = resolve_spin_num(input_file_name, 24);

    std::cout << std::setprecision(15);

    space::tools::Timer timer;
    timer.start();
    if (spin_num == 0) {
        space::multiThread::auto_multi_thread(single_thread_job<0>, work_dir, input_file, start_task_id, end_task_id,
                                              output_file, log_file);
    } else if (spin_num == 1) {
        space::multiThread::auto_multi_thread(single_thread_job<1>, work_dir, input_file, start_task_id, end_task_id,
                                              output_file, log_file);
    } else if (spin_num == 2) {
        space::multiThread::auto_multi_thread(single_thread_job<2>, work_dir, input_file, start_task_id, end_task_id,
                                              output_file, log_file);
    } else if (spin_num == 3) {
        space::multiThread::auto_multi_thread(single_thread_job<3>, work_dir, input_file, start_task_id, end_task_id,
                                              output_file, log_file);
    } else {
        std::cout << "wrong spin number!\n";
    }
    std::cout << "\r\n Time:" << timer.get_time() << '\n';
    return 0;
}*/

void get_line(std::fstream& is, std::string& str) {
    std::getline(is, str);
    if(!is)
      throw secular::StopFlag::eof;
}

void unpack_args(std::string  const& str, std::vector<double>& vec, bool DA, size_t spin_num) {
    std::stringstream is{str};

    size_t token_num = 22 + static_cast<size_t>(!DA) + spin_num*3;
    vec.reserve(token_num);

    double tmp;

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

    switch(token_num) {
        case 22 :
            return std::make_tuple(double_average, 0);
        case 23 :
            return std::make_tuple(single_average, 0);
        case 25 :
            return std::make_tuple(double_average, 1);
        case 26 :
            return std::make_tuple(single_average, 1);
        case 28 :
            return std::make_tuple(double_average, 2);
        case 29 :
            return std::make_tuple(single_average, 2);
        case 31 :
            return std::make_tuple(double_average, 3);
        case 32 :
            return std::make_tuple(single_average, 3);
        default :
            std::cout << "wrong input format!\n";
            throw secular::StopFlag::input_err;
    }
}

void single_thread_job(std::string work_dir, ConcurrentFile input, size_t start_id, size_t end_id, ConcurrentFile output) {
    //using Container = std::array<double, 12>;
    //using stepper_type = boost::numeric::odeint::runge_kutta_fehlberg78<Container>;
    double ini_dt = 0.1 * secular::consts::year;
    //secular::Task<spin_num> task;
    std::string entry;
    for (;;) {
        try{
            input.execute(get_line, entry);

            auto [DA, spin_num] = resolve_sim_type(entry);

            std::vector<double> v;

            unpack_args(entry, v, DA, spin_num);

            if(spin_num == 0) {

            } else if(spin_num == 1){

            } else if(spin_num == 2){

            } else if(spin_num == 3){

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

    //auto log_file = make_thread_safe_fstream(work_dir + "log.txt", std::fstream::out);

    std::cout << std::setprecision(15);

    space::tools::Timer timer;
    timer.start();
    space::multiThread::auto_multi_thread(single_thread_job, work_dir, input_file, start_task_id, end_task_id, output_file);
    std::cout << "\r\n Time:" << timer.get_time() << '\n';
    return 0;
}
