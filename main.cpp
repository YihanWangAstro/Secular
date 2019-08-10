
#include <iostream>
//#include "secular.h"
#include "refactor.h"
#include "SpaceHub/src/multi-thread/multi-thread.hpp"
#include "SpaceHub/src/tools/config-reader.hpp"
#include "SpaceHub/src/tools/timer.hpp"
#include "boost/numeric/odeint.hpp"
#include <iomanip>
#include <cstdlib>
#include <algorithm>
using namespace space::multiThread;

double int_error = 1e-13;

size_t resolve_spin_num(std::string const& file_name){
    std::fstream init(file_name, std::fstream::in);
    std::string line;
    std::getline(init, line);

    bool in_space = true;
    size_t token_num = 0;

    for(auto c : line){
        if(std::isspace(c)){
            in_space = true;
        } else if(in_space) {
            token_num ++;
            in_space = false;
        }
    }

    if(token_num == 22){
        return 0;
    } else if(token_num == 25) {
        return 1;
    } else if(token_num == 28) {
        return 2;
    } else if(token_num == 31) {
        return 3;
    } else {
        return 999;
    }
}

/*
template<size_t spin_num>
void single_thread_job(std::string work_dir, ConcurrentFile input, size_t start_id, size_t end_id, ConcurrentFile output, ConcurrentFile log) {
    using Sys = secular::Secular<spin_num, secular::Controler>;
    using stepper_type = boost::numeric::odeint::runge_kutta_fehlberg78<typename Sys::Container>;
    double ini_dt = 0.1 * secular::year;
    secular::Task<spin_num> task;
    for(;;) {
        if(input >> task) {
            if(start_id <= task.ctrl.id && task.ctrl.id <= end_id) {
                secular::deg_to_rad(task.obt_args);
                Sys eqns{work_dir, output, task.obt_args, task.ctrl};
                try{
                    boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::make_controlled(int_error ,int_error , stepper_type() ), eqns, eqns.initial_conds, 0.0, eqns.ctrl.end_time, ini_dt, eqns);
                } catch (secular::StopFlag flag) {
                    if(flag == secular::StopFlag::shrink) {
                    //
                    }
                }
            }
        } else {
            break;
        } 
  }
}
*/

template<typename Container, typename Obt>
void to_init(Container& c, Obt const& orbit) {
        Vec3d L1 = secular::calc_angular_mom(orbit.m1, orbit.m2, orbit.a_in) * sqrt(1 - orbit.e_in * orbit.e_in) * secular::unit_j(orbit.i_in, orbit.Omega_in);
        Vec3d L2 = secular::calc_angular_mom(orbit.m1 + orbit.m2, orbit.m3, orbit.a_out) * sqrt(1 - orbit.e_out * orbit.e_out) * secular::unit_j(orbit.i_out, orbit.Omega_out);
        Vec3d e1 = orbit.e_in * secular::unit_e(orbit.i_in, orbit.omega_in, orbit.Omega_in);
        Vec3d e2 = orbit.e_out * secular::unit_e(orbit.i_out, orbit.omega_out, orbit.Omega_out);
        c[0] = L1.x, c[1] = L1.y, c[2] = L1.z;
        c[3] = e1.x, c[4] = e1.y, c[5] = e1.z;
        c[6] = L2.x, c[7] = L2.y, c[8] = L2.z;
        c[9] = e2.x, c[10] = e2.y, c[11] = e2.z;
}

template<size_t spin_num>
void single_thread_job(std::string work_dir, ConcurrentFile input, size_t start_id, size_t end_id, ConcurrentFile output, ConcurrentFile log) {
    using Container = std::array<double,12>;
   
    using stepper_type = boost::numeric::odeint::runge_kutta_fehlberg78<Container>;
    double ini_dt = 0.1 * secular::year;
    secular::Task<spin_num> task;
    for(;;) {
        if(input >> task) {
            if(start_id <= task.ctrl.id && task.ctrl.id <= end_id) {
                secular::deg_to_rad(task.obt_args);

                std::array<double,12> inits;

                to_init(inits, task.obt_args);

                secular::SecularArg args{task.obt_args.m1, task.obt_args.m2, task.obt_args.m3};

                try{
                    //boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::make_controlled(int_error ,int_error , stepper_type() ), eqns, inits, 0.0, task.ctrl.end_time, ini_dt, eqns);
                    boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::make_controlled(int_error ,int_error , stepper_type() ), secular::serilize(args, secular::quad_kozai<secular::SecularArg, Container>), inits, 0.0, task.ctrl.end_time, ini_dt);
                } catch (secular::StopFlag flag) {
                    if(flag == secular::StopFlag::shrink) {
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

    //space::tools::ConfigReader config("config.txt");
    //int_error = config.get<double>("error");

    //reading configure file name, and input file name from command line
    std::string input_file_name;
    std::string work_dir;
    size_t start_task_id, end_task_id;

    space::tools::read_command_line(argc, argv, input_file_name, start_task_id, end_task_id, work_dir);

    const int dir_err = system( ("mkdir -p " + work_dir).c_str() );
    if (dir_err == -1)
    {
        std::cout << "Error creating directory!n";
        exit(1);
    }

    work_dir += "/";

    auto input_file = make_thread_safe_fstream(input_file_name, std::fstream::in);

    auto output_file = make_thread_safe_fstream(work_dir + "statistics.txt", std::fstream::out);

    auto log_file = make_thread_safe_fstream(work_dir + "log.txt", std::fstream::out);

    size_t spin_num = resolve_spin_num(input_file_name);

    std::cout << std::setprecision(15);

    space::tools::Timer timer;
    timer.start();
    if(spin_num == 0) {
        space::multiThread::auto_multi_thread(single_thread_job<0>, work_dir, input_file, start_task_id, end_task_id, output_file, log_file);
    } else if(spin_num == 1) {
        space::multiThread::auto_multi_thread(single_thread_job<1>, work_dir, input_file, start_task_id, end_task_id, output_file, log_file);
    } else if(spin_num == 2) {
        space::multiThread::auto_multi_thread(single_thread_job<2>, work_dir, input_file, start_task_id, end_task_id, output_file, log_file);
    } else if(spin_num == 3) {
        space::multiThread::auto_multi_thread(single_thread_job<3>, work_dir, input_file, start_task_id, end_task_id, output_file, log_file);
    } else {
        std::cout << "wrong spin number!\n";
    }
    std::cout << "\r\n Time:" << timer.get_time() << '\n';
    return 0;
}