#include "occ_context_probe.hpp"
#include <iostream>
#include <exception>
int main(int argc, char **argv) {
  if (argc != 2) return 2;
  try {
    run_occ_context_probe(argv[1]);
    std::cout << "Native OCC identity, reader-state, located/shared/coincident topology, mesh and lifetime checks passed.\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << "\n";
    return 1;
  }
}
