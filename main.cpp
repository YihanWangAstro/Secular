
#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include "SpaceHub/src/multi-thread/multi-thread.hpp"
#include "SpaceHub/src/tools/config-reader.hpp"
#include "SpaceHub/src/tools/timer.hpp"
#include "boost/numeric/odeint.hpp"
#include "observer.h"
#include "secular.h"

using namespace space::multi_thread;
using namespace secular;

double ATOL = 1e-13;
double RTOL = 1e-13;

bool get_line(std::fstream &is, std::string &str) {
  std::getline(is, str);
  if (!is)
    return false;
  else
    return true;
}

constexpr size_t ARGS_OFFSET = 3;
constexpr size_t PARAMETER_NUM = 25;

auto call_ode_int(std::string work_dir, ConcurrentFile output, secular::Controller const &ctrl,
                  std::vector<double> const &init_args) {
  using namespace boost::numeric::odeint;

  using Container = secular::SecularArray;

  auto [task_id, t_end, out_dt] =
      secular::cast_unpack<decltype(init_args.begin()), size_t, double, double>(init_args.begin());

  // space::display(std::cout, task_id, t_end, out_dt);

  auto const [m1, m2, m3, a_in_init] = secular::unpack_args<4>(init_args.begin() + ARGS_OFFSET);

  std::fstream f_out;

  if (secular::is_on(out_dt)) {
    f_out.open(work_dir + "secular_" + std::to_string(task_id) + ".txt", std::fstream::out);
    f_out << std::setprecision(12);
  }

  Container data;

  initialize_orbit_args(ctrl.ave_method, data, init_args.begin() + ARGS_OFFSET);

  secular::SecularConst const_parameters{m1, m2, m3};

  double dt = 0.1 * secular::consts::year;

  double time = 0;

  // auto stepper = make_controlled(ATOL, RTOL, runge_kutta_fehlberg78<Container>());

  auto stepper = bulirsch_stoer<Container>{ATOL, RTOL};

  secular::Stream_observer writer{f_out, out_dt};

  secular::SMA_Determinator stop{const_parameters.a_in_coef(), a_in_init * ctrl.GW_in_ratio};

  auto func = secular::Dynamic_dispatch<Container>(ctrl, const_parameters);
  writer(data, time);

  // STATIC_DISPATH(ctrl, const_parameters,

  for (; time <= t_end && !stop(data, time);) {
    constexpr size_t max_attempts = 500;

    controlled_step_result res = success;
    size_t trials = 0;
    do {
      res = stepper.try_step(func, data, time, dt);
      trials++;
    } while ((res == fail) && (trials < max_attempts));

    if (trials == max_attempts) {
      return ReturnFlag::max_iter;
    }
    writer(data, time);
  }
  //)

  output << PACK(task_id, ' ', time, ' ', data, "\r\n");
  output.flush();

  return ReturnFlag::finish;
}

void single_thread_job(Controller const &ctrl, std::string work_dir, ConcurrentFile input, ConcurrentFile output,
                       ConcurrentFile log) {
  std::string entry;
  for (; input.execute(get_line, entry);) {
    size_t task_id = static_cast<size_t>(std::stoi(entry));

    std::vector<double> v;

    secular::unpack_args_from_str(entry, v, PARAMETER_NUM);

    ReturnFlag res = call_ode_int(work_dir, output, ctrl, v);

    if (res == ReturnFlag::max_iter) {
      log << std::to_string(task_id) + ":Max iteration number reaches!\n";
      log.flush();
    }
  }
}

size_t decide_thread_num(std::string const &user_specified_core_num, std::string const &input_file_path) {
  std::ifstream input_file{input_file_path};

  size_t task_num = std::count(std::istreambuf_iterator<char>(input_file), std::istreambuf_iterator<char>(), '\n');

  size_t cpu_num = space::multi_thread::machine_thread_num;

  if (user_specified_core_num != "auto") {
    if (secular::is_number(user_specified_core_num)) {
      cpu_num = std::stoi(user_specified_core_num);
    } else {
      std::cout << "wrong format of the first argument(cpu core number)!\n";
      exit(0);
    }
  }

  return std::min(task_num, cpu_num);
}

int main(int argc, char **argv) {
  std::cout << logo << std::endl;
  std::ios::sync_with_stdio(false);
  std::string input_file_name;
  std::string cfg_file_name;
  std::string work_dir;
  std::string user_specified_core_num;

  space::tools::read_command_line(argc, argv, cfg_file_name);

  space::tools::ConfigReader cfg{cfg_file_name};

  secular::Controller ctrl{cfg};

  ATOL = cfg.get<double>("absolute_tolerance");

  RTOL = cfg.get<double>("relative_tolerance");

  work_dir = cfg.get<std::string>("output_dir");

  input_file_name = cfg.get<std::string>("input");

  user_specified_core_num = cfg.get<std::string>("cpu_num");

  size_t thread_num = decide_thread_num(user_specified_core_num, input_file_name);

  std::cout << thread_num << " thread(s) will be created for calculation." << std::endl;

  const int dir_err = system(("mkdir -p " + work_dir).c_str());
  if (dir_err == -1) {
    std::cout << "Error creating directory!\n";
    return 0;
  }

  work_dir += "/";

  auto input_file = make_thread_safe_fstream(input_file_name, std::fstream::in);

  auto output_file = make_thread_safe_fstream(work_dir + "last_state.txt", std::fstream::out);

  auto log_file = make_thread_safe_fstream(work_dir + "log.txt", std::fstream::out);

  log_file << secular::get_log_title(ctrl) + "\r\n";
  log_file.flush();

  space::tools::Timer timer;
  timer.start();
  space::multi_thread::multi_thread(thread_num, single_thread_job, ctrl, work_dir, input_file, output_file, log_file);
  std::cout << "\r\n Time:" << timer.get_time() << " s\n";
  return 0;
}
