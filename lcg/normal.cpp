#include <algorithm>
#include <math.h>
#include <iostream>

#include "normal.h"

// Inverse Function Method
double Generator::Normal::GetSpecial() {
  // x ~ N(0, 1)
  double x = 0;
  for (int i = 0; i < 12; ++i) {
    // Get eps ~ U[0, 1]
    x += Get();
  }
  x -= 6;
  // Transform x ~ N(0, 1) into x ~ N(mu, sigma2)
  x *= sigma;
  x += mu;
  return x;
}



int Generator::Normal::GetNumParametrs() {
  return 2;
}



double Generator::Normal::f(const double x) {
  double pow = -x * x / 2;
  return exp(pow) / sqrt(2.0 * acos(-1));
}



void Generator::Normal::GetTheorFrecuency(std::vector<double>& theor_num,
  const double min_val) {
  std::cout << "Normal distribution:" << std::endl;
  double cur_val = min_val;
  for (int i = 0; i < theor_num.size(); ++i) {
    double mid_val = (2 * cur_val + h) / 2;
    mid_val -= mu;
    mid_val /= sigma;
    theor_num[i] = h * num_elements * f(mid_val) / sigma;
    cur_val += h;
  }
}
