
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

static size_t cpu_core{0};

static int SA_flag{0};

static int oct_flag{0};

static int gr_flag{0};

static int gw_flag{0};

static int ll_flag{0};

static int SiLi_flag{0};

static int SiLo_flag{0};

static int SoLi_flag{0};

static int SoLo_flag{0};

static bool stop_flag{false};

std::string prefix;

static double stop_a{0.0};

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

int read_options(int argc, char** argv) {
  cpu_core == space::multi_thread::machine_thread_num;

  int opt_name;

  struct option long_options[] = {/* These options set a flag. */
                                  {"oct", no_argument, &oct_flag, 1},
                                  {"LiLo", no_argument, &ll_flag, 1},
                                  {"SiLi", no_argument, &SiLi_flag, 1},
                                  {"SiLo", no_argument, &SiLo_flag, 1},
                                  {"SoLi", no_argument, &SoLi_flag, 1},
                                  {"SoLo", no_argument, &SoLo_flag, 1},
                                  /* These options don’t set a flag. We distinguish them by their indices. */
                                  {"method", required_argument, 0, 'm'},
                                  {"gr", required_argument, 0, 'g'},
                                  {"gw", required_argument, 0, 'w'},
                                  {"o", required_argument, 0, 'o'},
                                  {"stop", required_argument, 0, 's'},
                                  {"np", required_argument, 0, 'n'},
                                  {0, 0, 0, 0}};

  while ((opt_name = getopt_long(argc, argv, "m:g:w:o:s:n:", long_options, 0)) != EOF) {
    switch (opt_name) {
      case 0:
        break;
      case 'm':
        set_optarg(SA_flag, optarg, {"da", "sa"}, std::array{0, 1});
        break;
      case 'g':
        set_optarg(gr_flag, optarg, {"io", "oi", "i", "o"}, std::array{3, 3, 1, 2});
        break;
      case 'w':
        set_optarg(gw_flag, optarg, {"io", "oi", "i", "o"}, std::array{3, 3, 1, 2});
        break;
      case 'o':
        prefix = std::string(optarg);
        break;
      case 's':
        stop_a = std::stod(optarg);
        stop_flag = true;
        break;
      case 'n':
        cpu_core = std::stoul(optarg);
        break;
      case '?':
        /* getopt_long already printed an error message. */
        break;
      default:
        abort();
    }
  }
  return option_idx;
}

bool get_line(std::fstream& is, std::string& str) {
  std::getline(is, str);
  if (!is)
    return false;
  else
    return true;
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

void run(ConcurrentFile input_file) {
  std::string entry;
  for (; input_file.execute(get_line, entry);) {
    std::vector<double> inits;
    str_to_vector(entry, inits);
  }
}

int main(int argc, char** argv) {
  read_options(argc, argv);

  std::string input_file_name{argv[optind]};

  auto input_file = make_thread_safe_fstream(input_file_name, std::fstream::in);

  return 0;
}
