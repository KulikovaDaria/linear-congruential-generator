#include <iostream>
#include <algorithm>
#include "uniform.h"

// Inverse Function Method
double Generator::Uniform::GetSpecial() {
  return a + (b - a) * Get();
}



int Generator::Uniform::GetNumParametrs() {
  return 2;
}



void Generator::Uniform::GetTheorFrecuency(std::vector<double>& theor_num,
  double min_val) {
  std::cout << "Uniform distribution:" << std::endl;
  fill(theor_num.begin(), theor_num.end(),
    num_elements * 1.0 / theor_num.size());
}
