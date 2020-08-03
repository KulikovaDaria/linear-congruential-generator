#include <map>
#include "special_lcg.h"

bool Generator::SpecialLCG::Test() {
  // Get a sequence val, with special distribution
  std::vector<double> val(num_elements, 0);
  for (int i = 0; i < num_elements; ++i) {
    val[i] = GetSpecial();
  }
  std::sort(val.begin(), val.end());
  if (is_discrete) {
    num_intervals = val[num_elements - 1] - val[0] + 1;
  }
  else {
    h = (val[num_elements - 1] - val[0]) / num_intervals;
  }
  // Devide into n intervals and calculate the frequencies num[i] for each one
  std::vector<long long> num(num_intervals, 0);
  for (int i = 0; i < num_elements; ++i) {
    int idx = (val[i] - val[0]) / h;
    idx = std::min(std::max(idx, 0), num_intervals - 1);
    ++num[idx];
  }
  // Theoretical frecuency
  std::vector<double> theor_num(num_intervals, 0);
  GetTheorFrecuency(theor_num, val[0]);
  return X2(num, theor_num, std::max(num_intervals - GetNumParametrs() - 1, 1));
}