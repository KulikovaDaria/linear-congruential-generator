#ifndef UNIFORM_H
#define UNIFORM_H

#include "special_lcg.h"

namespace Generator {
  class Uniform : public SpecialLCG {
  public:
    Uniform(const double a, const double b)
      :a(a), b(b) {
      is_discrete = 0;
    }
    virtual double GetSpecial() override;

  private:
    // Inherited via SpecialLCG
    virtual int GetNumParametrs() override;
    virtual void GetTheorFrecuency(std::vector<double>& theor_num,
      double min_val) override;

    const double a{0};
    const double b{1};
  };
}

#endif