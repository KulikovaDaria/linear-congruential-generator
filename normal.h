#ifndef NORMAL_H
#define NORMAL_H

#include "special_lcg.h"

namespace Generator {
  class Normal : public SpecialLCG {
  public:
    Normal(const double mu, const double sigma2)
      :mu(mu), sigma2(sigma2), sigma(sqrt(sigma2)) {
      is_discrete = 0;
      num_intervals = num_elements / 100;
    }
    virtual double GetSpecial() override;


  private:
    // Inherited via SpecialLCG
    virtual int GetNumParametrs() override;
    virtual void GetTheorFrecuency(std::vector<double>& theor_num,
      double min_val) override;
    double f(const double x);

    const double mu{0};
    const double sigma2{1};
    const double sigma{1};
  };
}

#endif
