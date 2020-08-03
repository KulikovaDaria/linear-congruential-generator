#ifndef BERNOULLI_H
#define BERNOULLI_H
#include "special_lcg.h"

namespace Generator {
  class Bernoulli : public SpecialLCG {
  public:
    Bernoulli(const double p)
      :p(p) {
      num_intervals = 2;
    }
    virtual double GetSpecial() override;

  private:
    // Inherited via SpecialLCG
    virtual int GetNumParametrs() override;
    virtual void GetTheorFrecuency(std::vector<double>& theor_num,
      double min_val) override;

    const double p{0.5};
  };
}

#endif