#include "geometric.h"

double Generator::Geometric::GetSpecial() {
  int x = 0;
  double eps = Get();
  double factor = p;
  double sum = p;
  while (sum < eps) {
    factor *= (1 - p);
    sum += factor;
    ++x;
  }
  return x;
}



int Generator::Geometric::GetNumParametrs() {
  return 1;
}



void Generator::Geometric::GetTheorFrecuency(std::vector<double>& theor_num,
  double min_val) {
  std::cout << "Geometric distribution:" << std::endl;
  double factor = p;
  for (int i = 0; i < num_intervals; ++i) {
    theor_num[i] = factor * num_elements;
    factor *= (1 - p);
  }
}
