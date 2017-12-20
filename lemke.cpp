#include <iostream>

#include "LCP.h"

int main() {

  // lemke's algorithm
  LCP problem;
  
  problem.read_LCP();
  
  problem.print_LCP();

  problem.lemke_algorithm();

  return 0;
}
