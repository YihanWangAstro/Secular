
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
    }
}

template<size_t spin_num>
void single_thread_job(std::string work_dir, ConcurrentFile input,  ConcurrentFile output, ConcurrentFile log) {
    using Sys = secular::Secular<spin_num, secular::Controler>;
    using stepper_type = boost::numeric::odeint::runge_kutta_fehlberg78<typename Sys::Container>;
    double ini_dt = 0.1 * secular::year;
    secular::Task<spin_num> task;
    for(;;) {
        if(input >> task) {
            //std::cout << task << '\n';
            secular::deg_to_rad(task.obt_args);
            Sys eqns{work_dir, output, task.obt_args, task.ctrl};
            try{
                boost::numeric::odeint::integrate_adaptive(boost::numeric::odeint::make_controlled(int_error ,int_error , stepper_type() ), eqns, eqns.initial_conds, 0.0, eqns.ctrl.end_time, ini_dt, eqns);
            } catch (secular::StopFlag flag) {
                if(flag == secular::StopFlag::shrink) {
                    //
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
    space::tools::read_command_line(argc, argv, input_file_name, work_dir);

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
        space::multiThread::auto_multi_thread(single_thread_job<0>, work_dir, input_file, output_file, log_file);
    } else if(spin_num == 1) {
        space::multiThread::auto_multi_thread(single_thread_job<1>, work_dir, input_file, output_file, log_file);
    } else if(spin_num == 2) {
        space::multiThread::auto_multi_thread(single_thread_job<2>, work_dir, input_file, output_file, log_file);
    } else if(spin_num == 3) {
        space::multiThread::auto_multi_thread(single_thread_job<3>, work_dir, input_file, output_file, log_file);
    } else {
        std::cout << "wrong spin number!\n";
    }
    std::cout << "\r\n Time:" << timer.get_time() << '\n';
    return 0;
}