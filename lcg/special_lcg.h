#ifndef SPECIAL_LCG_H
#define SPECIAL_LCG_H

#include <algorithm>
#include <iostream>
#include <math.h>
#include <vector>
#include "lcg.h"

namespace Generator {
  class SpecialLCG : protected LCG {
  public:
    bool Test();
  protected:
    virtual double GetSpecial() = 0;
    virtual int GetNumParametrs() = 0;
    virtual void GetTheorFrecuency(std::vector<double>& theor_num,
      double min_val) = 0;

    bool is_discrete{1};
    const int num_elements{100000};
    int num_intervals{int(1.5 + 3.322 * log(num_elements))};
    double h{1};
  };
}

#endif