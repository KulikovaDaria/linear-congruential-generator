#include <iomanip>
#include <iostream>
#include <vector>
#include <omp.h>

#include "lcg.h"
#include "normal.h"
#include "uniform.h"
#include "exponential.h"
#include "binomial.h"
#include "bernoulli.h"
#include "geometric.h"
#include "poisson.h"
#include "hypergeometric.h"



int main() {
  using namespace Generator;
  std::cout << std::setprecision(10) << std::fixed;

  std::vector<int> q{2, 3, 7, 11, 31, 151, 331};
  LCG g;

  std::cout << "a % p = " << g.a % g.m << std::endl;
  std::cout << " a ^ ((p-1)/q) % p = " << std::endl;
  for (int i = 0; i < q.size(); ++i) {
    std::cout << "q =  " << q[i] << "  =>  " << g.Pow(g.a, (g.m - 1) / q[i]) << std::endl;
  }
  std::cout << std::endl;
  g.Test();


 // long long s = g.a * g.a;
 // std::cout << s % g.m << std::endl;

 

/*  Bernoulli Bern(0.9);
  Bern.Test();
  Binomial Bin(5, 0.5);
  Bin.Test();
  Geometric Geom(0.1);
  Geom.Test();
  Hypergeometric HGeom(12, 10, 5);
  HGeom.Test();
  Poisson Pois(0.2);
  Pois.Test();

  Uniform U(-10.1, 8.6);
  U.Test();
  Normal N(-3, 2.2);
  N.Test();
  Exponential Exp(0.2);
  Exp.Test();*/

  return 0;
 }
