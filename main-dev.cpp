
#include <getopt.h>
#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include "SpaceHub/src/multi-thread/multi-thread.hpp"
#include "SpaceHub/src/tools/config-reader.hpp"
#include "SpaceHub/src/tools/timer.hpp"
#include "boost/numeric/odeint.hpp"
#include "observer.h"
#include "secular.h"

using namespace space::multi_thread;
using namespace secular;

double INT_ERROR = 1e-13;

static bool COLLECT_FLAG{true};

static bool STAT_FINAL_FLAG{false};

static bool STEP_OUTPUT_FLAG{true};

static size_t CPU_CORE{1};

static bool SINGLE_AVE_FLAG{false};

static bool OCT_FLAG{false};

static OrbIdx GR_FLAG{OrbIdx::off};

static OrbIdx GW_FLAG{OrbIdx::off};

static SLstat SL_STAT{deS::off, deS::off, deS::off, deS::off, deS::off};

static bool STOP_FLAG{false};

std::string prefix{"./"};

static double STOP_A{0.0};

template <typename Container>
struct Equation_of_Motion {
  Equation_of_Motion(SecularConst const& _args) : args{&_args} {}

  void operator()(Container const& x, Container& dxdt, double t) {
    std::fill(dxdt.begin(), dxdt.end(), 0);

    Lidov_Kozai(SINGLE_AVE_FLAG, OCT_FLAG, *args, x, dxdt, t);

    GR_precession(SINGLE_AVE_FLAG, GR_FLAG, *args, x, dxdt, t);

    GW_radiation(GW_FLAG, *args, x, dxdt, t);

    deSitter_precession(SINGLE_AVE_FLAG, SL_STAT, *args, x, dxdt, t);
  }
  SecularConst const* args;
};

/*
size_t get_spin_num() {
  size_t spin_num = 0;
  if (SiLi_flag == 1 || SiLo_flag == 1) {
    spin_num += 2;
  }
  if (SoLi_flag == 1 || SoLo_flag == 1) {
    spin_num += 1;
  }
}*/

template <typename Arg, size_t sz>
void set_optarg(Arg& val, char const* optarg, std::array<char const*, sz> const& keys,
                std::array<Arg, sz> const& values) {
  size_t i = 0;
  for (auto key : keys) {
    if (strcmp(key, optarg) == 0) {
      val = values[i];
      return;
    }
    i++;
  }

  printf("Undefined option value '%s', use", optarg);
  for (auto key : keys) {
    printf(" '%s'", key);
  }
  printf("instead.\n");

  exit(0);
}

static const std::string help_info{
    "Usage: secular [options] input_file\n\
    Options : (the first possible value is the default value.)\n\
    --oct       'off', 'on'                    Octupole Lidov-Kozai effect.\n\
    --method    'double', 'single'             Orbital average method.\n\
    --gr        'off', 'in', 'out', 'both'     GR precession of inner/outer orbit.\n\
    --gw        'off', 'in', 'out', 'both'     Gravitational radiation of inner/outer orbit.\n\
    --LiLo      'off', 'both', 'noback'        (L_in, L_out) coupling and feedback reaction.\n\
    --SiLi      'off', 'both', 'noback'        (S_in, L_in) coupling and feedback reaction.\n\
    --SiLo      'off', 'both', 'noback'        (S_in, L_out) coupling and feedback reaction.\n\
    --SoLi      'off', 'both', 'noback'        (S_out, L_in) coupling and feedback reaction.\n\
    --SoLo      'off', 'both', 'noback',       (S_out, L_out) coupling and feedback reaction.\n\
    --o          <string>                      The output path.\n\
    --np         <number>                      Number of threads. The default value is the number of the logical CPU cores.\n\
    --stop       <number>                      The termination semi-major axis of the inner orbit when any decay effect is turned on. [in length unit of input file]\n\
    --collect   'off', 'on'                    Collect the end status of each job into a file."};

void read_options(int argc, char** argv) {
  // cpu_core == space::multi_thread::machine_thread_num;
  int opt_name;

  struct option long_options[] = {/* These options set a flag. */
                                  {"help", no_argument, 0, 'h'},       {"oct", required_argument, 0, 'b'},
                                  {"LiLo", required_argument, 0, '1'}, {"SiLi", required_argument, 0, '2'},
                                  {"SiLo", required_argument, 0, '3'}, {"SoLi", required_argument, 0, '4'},
                                  {"SoLo", required_argument, 0, '5'}, {"method", required_argument, 0, 'm'},
                                  {"gr", required_argument, 0, 'g'},   {"gw", required_argument, 0, 'w'},
                                  {"o", required_argument, 0, 'o'},    {"stop", required_argument, 0, 's'},
                                  {"np", required_argument, 0, 'n'},   {"collect", required_argument, 0, 'c'},
                                  {"tol", required_argument, 0, 't'},  {0, 0, 0, 0}};

  while ((opt_name = getopt_long(argc, argv, "m:g:w:o:s:n:1:2:3:4:5:b:", long_options, 0)) != EOF) {
    switch (opt_name) {
      case 0:
        break;
      case 'h':
        std::cout << help_info << std::endl;
        exit(0);
        break;
      case 'm':
        set_optarg(SINGLE_AVE_FLAG, optarg, {"double", "single"}, std::array{false, true});
        break;
      case 'b':
        set_optarg(OCT_FLAG, optarg, {"on", "off"}, std::array{true, false});
        break;
      case 'c':
        set_optarg(COLLECT_FLAG, optarg, {"on", "off"}, std::array{true, false});
        break;
      case 'g':
        set_optarg(GR_FLAG, optarg, {"both", "in", "out", "off"},
                   std::array{OrbIdx::in_out, OrbIdx::in, OrbIdx::out, OrbIdx::off});
        break;
      case 'w':
        set_optarg(GW_FLAG, optarg, {"both", "in", "out", "off"},
                   std::array{OrbIdx::in_out, OrbIdx::in, OrbIdx::out, OrbIdx::off});
        break;
      case '1':
        set_optarg(SL_STAT.LL, optarg, {"both", "noback", "off"}, std::array{deS::bc, deS::on, deS::off});
        break;
      case '2':
        set_optarg(SL_STAT.Sin_Lin, optarg, {"both", "noback", "off"}, std::array{deS::bc, deS::on, deS::off});
        break;
      case '3':
        set_optarg(SL_STAT.Sin_Lout, optarg, {"both", "noback", "off"}, std::array{deS::bc, deS::on, deS::off});
        break;
      case '4':
        set_optarg(SL_STAT.Sout_Lin, optarg, {"both", "noback", "off"}, std::array{deS::bc, deS::on, deS::off});
        break;
      case '5':
        set_optarg(SL_STAT.Sout_Lout, optarg, {"both", "noback", "off"}, std::array{deS::bc, deS::on, deS::off});
        break;
      case 'o':
        prefix = std::string(optarg);
        break;
      case 's':
        STOP_A = std::stod(optarg);
        STOP_FLAG = true;
        break;
      case 'n':
        CPU_CORE = static_cast<size_t>(std::stoi(optarg));
        break;
      case 't':
        INT_ERROR = pow(10.0, -std::stoi(optarg));
        break;
      case '?':
        /* getopt_long already printed an error message. */
        break;
      default:
        abort();
    }
  }
}

void str_to_vector(std::string const& str, std::vector<double>& vec) {
  std::stringstream is{str};
  vec.reserve(30);
  double tmp;
  for (; is;) {
    is >> tmp;
    vec.emplace_back(tmp);
  }
}

inline size_t proper_input_args_number() { return 18; }

void run(ConcurrentFile input_file, ConcurrentFile stat_output) {
  using namespace boost::numeric::odeint;

  using Container = secular::SecularArray;

  std::string entry;
  for (; input_file.execute(get_line, entry);) {
    std::vector<double> inits;

    str_to_vector(entry, inits);

    if (inits.size() < proper_input_args_number()) {
      std::cout << "wrong arguments number in line: " << entry << std::endl;
      continue;
    }
    size_t task_id = static_cast<size_t>(inits[0]);
    double t_end = inits[1];
    double out_dt = inits[2];
    double m1 = inits[3];
    double m2 = inits[4];
    double m3 = inits[5];

    size_t spin_num = (inits.size() - 18) / 3;

    Container data{spin_num};

    for (size_t i = 6; i < inits.size(); ++i) {
      data[i - 6] = inits[i];
    }

    secular::SecularConst const_parameters{m1, m2, m3};

    std::cout << task_id << ' ' << t_end << ' ' << out_dt << ' ' << m1 << ' ' << ' ' << m2 << ' ' << spin_num << ' '
              << data << std::endl;

    double dt = t_end / 10000;

    double time = 0;

    std::fstream f_out;

    if (STEP_OUTPUT_FLAG) {
      f_out.open(prefix + "output_" + std::to_string(task_id) + ".txt", std::fstream::out);
      f_out << std::setprecision(12);
    }

    secular::Stream_observer writer{f_out, out_dt, STEP_OUTPUT_FLAG};

    secular::SMA_Determinator stop{const_parameters.a_in_coef(), STOP_A, STOP_FLAG};

    auto func = Equation_of_Motion<Container>(const_parameters);

    writer(data, time);

    auto stepper = make_controlled(INT_ERROR, INT_ERROR, runge_kutta_fehlberg78<Container>());

    auto stepper = bulirsch_stoer<Container>{INT_ERROR, INT_ERROR};

    for (; time <= t_end && !stop(data, time);) {
      constexpr size_t max_attempts = 500;

      controlled_step_result res = success;
      size_t trials = 0;
      do {
        res = stepper.try_step(func, data, time, dt);
        trials++;
      } while ((res == fail) && (trials < max_attempts));

      if (trials == max_attempts) {
        std::cout << "Reach max iteration number in task: " << task_id << std::endl;
        continue;
      }
      writer(data, time);
    }

    stat_output << PACK(task_id, ' ', time, ' ', data, "\r\n");
    stat_output.flush();
  }
}

int main(int argc, char** argv) {
  read_options(argc, argv);

  std::string input_file_name{argv[optind]};

  std::cout << input_file_name << ' ' << CPU_CORE << std::endl;

  auto input_file = make_thread_safe_fstream(input_file_name, std::fstream::in);

  auto output_file = make_thread_safe_fstream(prefix + "/collection.txt", std::fstream::out);

  multi_thread(CPU_CORE, run, input_file, output_file);

  return 0;
}
