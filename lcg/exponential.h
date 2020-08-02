#ifndef EXPONENTIAL_H
#define EXPONENTIAL_H

#include "special_lcg.h"

namespace Generator {
  class Exponential : public SpecialLCG {
  public:
    Exponential(const double lambda)
      :lambda(lambda) {
      is_discrete = 0;
    }
    virtual double GetSpecial() override;

  private:
    // Inherited via SpecialLCG
    virtual int GetNumParametrs() override;
    virtual void GetTheorFrecuency(std::vector<double>& theor_num,
      double min_val) override;
    double f(const double x);

    const double lambda{1};
  };
}

#endif