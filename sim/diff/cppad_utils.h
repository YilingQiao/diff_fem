#ifndef CPPAD_UTILS_H
#define CPPAD_UTILS_H

#include <assert.h>
#include <stdio.h>
#include <string>


#include <Eigen/Core>
#include <cppad/cppad.hpp>
#include "math.h"

/**
 * Supporting functions for automatic differentiation with CppAD.
 * https://coin-or.github.io/CppAD/doc/cppad.htm
 */
template <typename InnerScalar = double>
struct CppADUtils {
  typedef typename CppAD::AD<InnerScalar> Scalar;

  static Scalar zero() { return Scalar(0.); }
  static Scalar one() { return Scalar(1.); }

  static Scalar two() { return Scalar(2.); }
  static Scalar half() { return Scalar(0.5); }
  static Scalar pi() { return Scalar(M_PI); }
  static Scalar half_pi() { return Scalar(M_PI / 2.); }

  static bool getBool(bool v) {  return v; }

  static Scalar scalar_from_double(double d) {
    return Scalar(d);
  }

  static Scalar scalar_from_string(const std::string &txt) {
    Scalar result = Scalar(atof(txt.c_str()));
    return result;
  }


  template <class T>
  static T cos1(const T& v) {
    using std::cos;
    return cos(v);
  }

  template <class T>
  static T sin1(const T& v) {
    using std::sin;
    return sin(v);
  }

  template <class T>
  static T abs1(const T& v) {
    using std::abs;
    return abs(v);
  }

  template <class T>
  static T exp1(const T& v) {
    using std::exp;
    return exp(v);
  }

  template <class T>
  static T sqrt1(const T& v) {
    using std::sqrt;
    return sqrt(v);
  }

  template <class T>
  static T atan2(const T& dy, const T& dx) {
    using std::atan2;
    return atan2(dy, dx);
  }

  static double getDouble(const Scalar& v) { return CppAD::Value(v); }
  // template <class T>
  // static double getDouble(const T& v) {
  //   return CppAD::Value(v);
  // }

  template <class T>
  static Scalar convert(T) = delete;  // C++11

  static Scalar convert(int value) { return Scalar(double(value)); }

  // template <class T>
  // static Scalar fraction(T, T) = delete;  // C++11

  static Scalar fraction(int num, int denom) {
    return Scalar(double(num) / double(denom));
  }

  static void FullAssert(bool a) {
    if (!a) {
      printf("!");
      assert(0);
      exit(0);
    }
  }
};

namespace std {
    template<>  CppAD::AD<double>  numeric_limits<CppAD::AD<double>>::min() {
      return CppAD::AD<double>(2.22507e-308);
    }
}

namespace Eigen {
  template<> struct NumTraits<CppAD::AD<double>> : GenericNumTraits<CppAD::AD<double>>
  {
    typedef CppAD::AD<double> Real;
    typedef CppAD::AD<double> NonInteger;
    typedef CppAD::AD<double> Nested;
 
    static inline Real epsilon() { return 2.22045e-16; }
    static inline Real min() { return 2.22507e-308; }
 
    enum {
      IsComplex = 0,
      IsInteger = 0,
      IsSigned = 1,
      RequireInitialization = 1,
      ReadCost = 1,
      AddCost = 3,
      MulCost = 3
    };
  };
 
}

#endif  // CPPAD_UTILS_H
