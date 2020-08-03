#ifndef POISSON_H
#define POISSON_H

#include "special_lcg.h"

namespace Generator {
  class Poisson : public SpecialLCG {
  public:
    Poisson(const double lambda)
      :lambda(lambda) {}
    virtual double GetSpecial() override;

  private:
    // Inherited via SpecialLCG
    virtual int GetNumParametrs() override;
    virtual void GetTheorFrecuency(std::vector<double>& theor_num,
      double min_val) override;

    const double lambda{1};
  };
}

#endif