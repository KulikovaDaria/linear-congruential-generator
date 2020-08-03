#include <math.h>
#include "binomial.h"

// Inverse Function Method
double Generator::Binomial::GetSpecial() {
  int x = 0;
  for (int i = 0; i < n; ++i) {
    if (Get() < p) {
      ++x;
    }
  }
  return x;
}



int Generator::Binomial::GetNumParametrs() {
  return 2;
}



int Generator::Binomial::C(const int n, const int k) {
  if (k > n - k) {
    return C(n, n - k);
  }
  int res = 1;
  for (int i = 1; i <= k; ++i) {
    res *= (n - k + i);
    res /= i;
  }
  return res;
}



double Generator::Binomial::f(const double x) {
  return C(n, x) * pow(p, x) * pow(1 - p, n - x);
}



void Generator::Binomial::GetTheorFrecuency(std::vector<double>& theor_num,
  double min_val) {
  std::cout << "Binomial distribution:" << std::endl;
  for (int i = 0; i <= n; ++i) {
    theor_num[i] = f(i) * num_elements;
  }
}
