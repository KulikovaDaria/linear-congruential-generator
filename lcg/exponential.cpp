#include "exponential.h"
#include <math.h>
#include "exponential.h"

// Inverse Function Method
double Generator::Exponential::GetSpecial() {
  return -log(Get()) / lambda;
}



int Generator::Exponential::GetNumParametrs() {
  return 1;
}



double Generator::Exponential::f(const double x) {
  return exp(-lambda * x);
}



void Generator::Exponential::GetTheorFrecuency(std::vector<double>& theor_num,
  double min_val) {
  std::cout << "Exponential distribution:" << std::endl;
  double left = min_val;
  double right = h;
  for (int i = 0; i < theor_num.size() - 1; ++i) {
    theor_num[i] = f(left) - f(right);
    theor_num[i] *= num_elements;
    left = right;
    right += h;
  }
  theor_num[theor_num.size() - 1] = f(left) * num_elements;
}
