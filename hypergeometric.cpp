#include "hypergeometric.h"

double Generator::Hypergeometric::GetSpecial() {
  int white = 0;
  int black = 0;
  double p = w  * 1.0 / (w + b);
  for (int i = 1; i <= n && white < w && black < b; ++i) {
    if (Get() < p) {
      ++white;
    }
    else {
      ++black;
    }
    p = (w - white) * 1.0 / (w + b - i);
  }
  if (white + black < n && white < w) {
    white = n - black;
  }
  return white;
}



int Generator::Hypergeometric::GetNumParametrs() {
  return 3;
}



long long Generator::Hypergeometric::C(const int n, const int k) {
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



double Generator::Hypergeometric::f(const int x) {
  if (n - x > b) {
    return 0;
  }
  return C(w, x) * C(b, n - x) * 1.0 / C(w + b, n);
}



void Generator::Hypergeometric::GetTheorFrecuency(std::vector<double>& theor_num,
  double min_val) {
  std::cout << "Hypergeometric distribution:" << std::endl;
  theor_num[0] = f(min_val) * num_elements;
  for (int i = 1; i < num_intervals; ++i) {
    theor_num[i] = theor_num[i - 1] * (w - min_val - i + 1)
      * (n - min_val - i + 1) / ((min_val + i) * (b - n + min_val + i));
  }
}
