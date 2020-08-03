cmake_minimum_required(VERSION 3.6)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER "cmake")

project( lcg )

add_executable ( lcg bernoulli.h bernoulli.cpp binomial.h binomial.cpp exponential.h exponential.cpp geometric.h geometric.cpp hypergeometric.h hypergeometric.cpp lcg.h lcg.cpp normal.h normal.cpp poisson.h poisson.cpp special_lcg.h special_lcg.cpp test.cpp uniform.h uniform.cpp )