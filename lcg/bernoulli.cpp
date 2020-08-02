#include "bernoulli.h"

// Inverse Function Method
double Generator::Bernoulli::GetSpecial() {
  return Get() < p;
}



int Generator::Bernoulli::GetNumParametrs() {
  return 1;
}



void Generator::Bernoulli::GetTheorFrecuency(std::vector<double>& theor_num,
  double min_val) {
  std::cout << "Bernoulli distribution:" << std::endl;
  theor_num[0] = (1 - p) * num_elements;
  theor_num[1] = p * num_elements;
}
