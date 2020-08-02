#include "poisson.h"

// Inverse Function Method
double Generator::Poisson::GetSpecial() {
  int x = 0;
  double mult = Get();
  double bound = exp(-lambda);
  while (mult >= bound) {
    ++x;
    mult *= Get();
  }
  return x;
}



int Generator::Poisson::GetNumParametrs() {
  return 1;
}



void Generator::Poisson::GetTheorFrecuency(std::vector<double>& theor_num,
  double min_val) {
  std::cout << "Poisson distribution:" << std::endl;
  theor_num[0] = exp(-lambda) * num_elements;
  for (int i = 1; i < num_intervals; ++i) {
    theor_num[i] = theor_num[i - 1] * lambda / i;
  }
}
