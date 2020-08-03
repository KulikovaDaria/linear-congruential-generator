#ifndef HYPERGEOMETRIC_H
#define HYPERGEOMETRIC_H

#include "special_lcg.h"

namespace Generator {
  class Hypergeometric : public SpecialLCG {
  public:
    Hypergeometric(const int w, const int b, const int n)
      :w(w), b(b), n(n) {}
    virtual double GetSpecial() override;

  private:
    // Inherited via SpecialLCG
    virtual int GetNumParametrs() override;
    virtual void GetTheorFrecuency(std::vector<double>& theor_num,
      double min_val) override;
    long long C(const int n, const int k);
    double f(const int x);

    const int w{1};
    const int b{1};
    const int n{1};
  };
}

#endif
