#ifndef BINOMIAL_H
#define BINOMIAL_H

#include "special_lcg.h"

namespace Generator {
  class Binomial : public SpecialLCG {
  public:
    Binomial() = default;
    Binomial(const int n, const double p)
      :n(n), p(p) {}
    virtual double GetSpecial() override;

  private:
    // Inherited via SpecialLCG
    virtual int GetNumParametrs() override;
    virtual void GetTheorFrecuency(std::vector<double>& theor_num,
      double min_val) override;
    int C(const int n, const int k);
    double f(const double x);

    const int n{1};
    const double p{1};
  };
}

#endif