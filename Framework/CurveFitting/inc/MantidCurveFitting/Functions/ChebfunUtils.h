#ifndef MANTID_CURVEFITTING_CHEBFUNUTILS_H_
#define MANTID_CURVEFITTING_CHEBFUNUTILS_H_

#include "MantidCurveFitting/Functions/ChebfunBase.h"
#include "MantidCurveFitting/GSLMatrix.h"
#include "MantidCurveFitting/GSLVector.h"
#include "MantidCurveFitting/Functions/SimpleChebfun.h"

namespace Mantid {
  namespace CurveFitting {
    namespace Functions {

      /// Create an interpolating matrix.
      GSLMatrix createInterpolatingMatrix(const ChebfunBase& base, std::vector<double> &x, bool isZeroOutside = true);
      GSLVector createBasisVector(const ChebfunBase& base, double x);
      ChebfunBase_sptr uniformFit(double startX, double endX, const std::vector<double> &y, std::vector<double> &p);
      ChebfunBase_sptr nonUniformFit(const std::vector<double> &x, const std::vector<double> &y, std::vector<double> &p, size_t m = 0);
      ChebfunBase_sptr inverse(const ChebfunBase& cheb, const std::vector<double> &y, std::vector<double> &p);
    }
  }
}

#endif // MANTID_CURVEFITTING_CHEBFUNUTILS_H_