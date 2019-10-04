
#include <iostream>
#include "SpaceHub/src/tools/config-reader.hpp"
#include "secular.h"

int main(int argc, char **argv) {
  std::ios::sync_with_stdio(false);

  std::string cfg_file_name;

  space::tools::read_command_line(argc, argv, cfg_file_name);

  space::tools::ConfigReader cfg{cfg_file_name};

  secular::Controller ctrl{cfg};

  std::cout << "\033[1;31mInput column names:\033[0m\n\033[32m" << ctrl.initial_format()
            << "\033[0m\n\n\n\033[1;31mOutput column names:\033[0m\n\033[32m" << ctrl.output_format() << "\033[0m"
            << std::endl;
  return 0;
}
