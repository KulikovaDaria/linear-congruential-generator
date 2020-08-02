#ifndef LCG_H
#define LCG_H

#include <vector>

namespace Generator {
  class LCG {
  public:
    LCG() = default;
    inline double Get();
    double GetTransition(const long long k);
    bool Compare(const long long k);
    bool Test();
    void ParallelTest();

//  protected:
    bool X2(const std::vector<long long>& v, const std::vector<double>& p,
      const long long degrees_of_freedom);

 // private:
    inline long long Add(const long long num1, const long long num2);
    inline long long Subtr(const long long num1, const long long num2);
    inline long long Mult(const long long num1, const long long num2);
    inline long long Div(const long long num1, const long long num2);
    inline long long Pow(const long long num, const long long pow);
    inline double GetEps(const long long x);
    inline long long GetX(const long long from);
    long long GetXTransition(const long long from, const long long k);
    double GoldsteinApproximation(const long long degrees_of_freedom);
    bool EquidistributionTest();
    bool EquidistributionTestParallel();
    bool SerialTest();
    bool SerialTestParallel();
    bool PokerTest();
    bool PokerTestParallel();
    bool CouponCollectorsTest();
    bool CouponCollectorsTestParallel();

    const long long n{1ll << 31};
    const long long m{(1ll << 31) - 1};
    long long a{62089917};
    const long long c{0};
    const long long x0{1};
    long long x{1};
  };
}

#endif
