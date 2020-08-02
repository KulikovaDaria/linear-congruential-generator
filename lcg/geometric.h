#ifndef GEOMETRIC_H
#define GEOMETRIC_H

#include "special_lcg.h"

namespace Generator {
  class Geometric : public SpecialLCG {
  public:
    Geometric(const double p)
      :p(p) {}
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