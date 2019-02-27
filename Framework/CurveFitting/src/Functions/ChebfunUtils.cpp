#include "MantidCurveFitting/Functions/ChebfunUtils.h"
#include "MantidCurveFitting/Functions/IntervalChebfun.h"

#include "MantidAPI/Algorithm.h"
#include "MantidAPI/WorkspaceProperty.h"
#include "MantidAPI/WorkspaceFactory.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidCurveFitting/Functions/CubicSpline.h"
#include "MantidCurveFitting/GSLVector.h"

#include <boost/lexical_cast.hpp>
#include <numeric>

#pragma warning(disable: 4996)

namespace Mantid {
namespace CurveFitting {
namespace Functions {

double *begin(GSLVector& vec) {
  return &vec[0];
}

double *end(GSLVector& vec) {
  return &vec[vec.size() - 1];
}

/// Vector-matrix multiplication: Tr(v) * M
GSLVector operator*(const GSLVector& v, const GSLMatrix& M) {
  if (v.size() != M.size1()) throw std::runtime_error("Cannot multiply vector by matrix: size mismatch.");
  GSLVector res(M.size2());
  for (size_t i = 0; i < res.size(); ++i) {
    double tmp = 0.0;
    for (size_t j = 0; j < v.size(); ++j) {
      tmp += v.get(j) * M.get(j, i);
    }
    res.set(i, tmp);
  }
  return res;
}

void setRow(GSLMatrix& M, size_t j, const GSLVector& row) {
  if (row.size() != M.size1()) throw std::runtime_error("Cannot set row: size mismatch.");
  for (size_t i = 0; i < row.size(); ++i) {
    M.set(i, j, row[j]);
  }
}

/// Create a matrix (rectangular in general) which if multiplied by a y-point vector
/// produces a vector of y-points calculated at given values of x.
/// @param x :: A vector of x-points where a function will be interpolated.
/// @return :: The interpolating matrix.
GSLMatrix createInterpolatingMatrix(const ChebfunBase& base, std::vector<double> &x, bool isZeroOutside)
{
  const size_t m = x.size();
  const size_t n = base.size();
  GSLMatrix M(m, n);
  auto& xp = base.xPoints();
  for (size_t i = 0; i < m; ++i)
  {
    const double xi = x[i];
    if (xi < base.startX() || xi > base.endX())
    {
      if (isZeroOutside)
      {
        for (size_t j = 0; j < n; ++j)
          M.set(i, j, 0.0);
      }
      else
        throw std::runtime_error("Cannot interpolate outside function domain.");
    }
    else {
      for (size_t j = 0; j < n; ++j)
      {
        const double xj = xp[j];
        double d = 1.0;
        for (size_t k = 0; k < n; ++k)
        {
          if (k == j) continue;
          const double xk = xp[k];
          d *= (xi - xk) / (xj - xk);
        }
        M.set(i, j, d);
      }
    }
  }
  return M;
}

/// Create a vector (the same size as the base) which if multiplied by a y-point vector
/// produces a function value at x.
/// @param base :: A chebfun base.
/// @param x :: An x-point where a function will be evaluated.
/// @return :: The basis matrix.
GSLVector createBasisVector(const ChebfunBase& base, double x) {
  GSLVector v(base.size());
  auto &xp = base.xPoints();
  if (x < base.startX() || x > base.endX()) {
    throw std::runtime_error("Cannot interpolate outside function domain.");
  }
  else {
    auto n = base.size();
    for (size_t j = 0; j < n; ++j) {
      const double xj = xp[j];
      double d = 1.0;
      for (size_t k = 0; k < n; ++k) {
        if (k == j)
          continue;
        const double xk = xp[k];
        d *= (x - xk) / (xj - xk);
      }
      v[j] = d;
    }
  }
  return v;
}

/**
  * Create a matrix that performs convolution when multiplied by a y-vector of another base.
  * c = M * y
  * @param r :: y-points in this base that define a convolution kernel.
  * @param base :: Another base. Y-points in that base can be multiplied by the created matrix
  * to calculate a convolution.
  * @return :: The convolution matrix.
  */
GSLMatrix createConvolutionMatrix(const ChebfunBase &kernel, const std::vector<double> &r, const ChebfunBase &base)
{
  if (r.size() != kernel.size()) throw std::invalid_argument("Cannot create convolution matrix: kernel function has wrong size.");

  // rw_j = r_j * w_j
  GSLVector rw = [&] {
    GSLVector w(kernel.integrationWeights());
    std::transform(begin(w), end(w), r.begin(), begin(w), std::multiplies<double>());
    return w;
  }();

  const size_t m = base.size();
  GSLMatrix M(m, m);

  for (size_t i = 0; i < m; ++i)
  {
    const double xi = base.xPoints()[i];
    std::vector<double> xj = kernel.xPoints();
    // xj = xi - xj
    std::transform(xj.begin(), xj.end(), xj.begin(), [&xi](double x) {return xi - x; });

    GSLMatrix S = createInterpolatingMatrix(base, xj, false);
    GSLVector row = rw * S;
    setRow(M, i, row);

  }

  return M;
}

ChebfunBase_sptr uniformFit(double startX, double endX, const std::vector<double> &y, std::vector<double> &p)
{
  auto n = y.size();
  ChebfunBase base(n - 1, -1, 1);

  // TODO: increase size of base1 iteratively until y[i] can be reproduces to required
  // accuracy. Accuracy can be an optional input argument.
  auto base1 = boost::make_shared<ChebfunBase>(4 * n - 1, startX, endX);
  auto &x = base1->xPoints();
  p.resize(x.size());
  const double D = endX - startX;
  for (size_t i = 0; i < x.size(); ++i)
  {
    double xx = x[i];
    double xx1 = -cos(M_PI*(xx - startX) / D);
    p[i] = base.eval(xx1, y);
  }
  return base1;
}

std::vector<double> roughRoots(const ChebfunBase& base, const std::vector<double> &p, double level)
{
  std::vector<double> rs;
  if (p.empty()) return rs;
  auto &x = base.xPoints();
  auto y1 = p.front() - level;
  for (size_t i = 1; i < p.size(); ++i) {
    auto y = p[i] - level;
    if (y == 0.0) {
      rs.push_back(x[i]);
    }
    else if (y1 * y < 0.0) {
      rs.push_back((-x[i - 1] * y + x[i] * y1) / (y1 - y));
    }
    y1 = y;
  }
  return rs;
}

ChebfunBase_sptr inverse(const ChebfunBase& cheb, const std::vector<double> &y, std::vector<double> &p)
{
  // TODO: increase accuracy. Use more accurate roots, iteratively increase n.
  // Check (approximately) that cheb is monotone
  bool increasing = y.back() - y.front() > 0.0;
  for (size_t i = 1; i < y.size(); ++i)
  {
    if (y[i] < y[i - 1])
    {
      if (increasing)
      {
        std::string msg = "Cannot find inverse: function is not increasing.\n" +
          boost::lexical_cast<std::string>(y[i - 1]) + ", " + boost::lexical_cast<std::string>(y[i]);
        throw std::runtime_error(msg);
      }
    }
    else if (!increasing) throw std::runtime_error("Cannot find inverse: function is not monotone.");
  }
  const double startY = increasing ? y.front() : y.back();
  const double endY = increasing ? y.back() : y.front();
  size_t n = cheb.order();
  auto base = boost::make_shared<ChebfunBase>(n, startY, endY);
  auto &yi = base->xPoints();
  p.resize(yi.size());
  p.front() = increasing ? cheb.startX() : cheb.endX();
  p.back() = increasing ? cheb.endX() : cheb.startX();
  for (size_t i = 1; i < yi.size() - 1; ++i)
  {
    auto roots = roughRoots(cheb, y, yi[i]);
    if (roots.size() != 1) throw std::runtime_error("Multiple or no roots.");
    p[i] = roots.front();
  }
  return base;
}

ChebfunBase_sptr nonUniformFit(const std::vector<double> &x, const std::vector<double> &y, std::vector<double> &p, size_t m)
{
  if (x.size() != y.size()) throw std::runtime_error("nonUniformFit: different vector sizes.");
  if (x.size() < 3) throw std::runtime_error("nonUniformFit: vectors must have at least 3 elements.");
  auto n = y.size() - 1;

  //ChebfunBase_sptr xInverseBase;
  //std::vector<double> xp;
  //{
  //  std::vector<double> xd;
  //  ChebfunBase_sptr xDirect = uniformFit(0, n, x, xd);
  //  xInverseBase = inverse(*xDirect, xd, xp);
  //}

  ChebfunBase base(n, -1, 1);
  if (m == 0) m = 2 * n;
  // TODO: increase size of base1 iteratively until y[i] can be reproduces to required
  // accuracy. Accuracy can be an optional input argument.
  auto fitBase = boost::make_shared<ChebfunBase>(m, x.front(), x.back());
  auto &b = fitBase->xPoints();
  p.resize(b.size());
  for (size_t i = 0; i < b.size(); ++i)
  {
    //double xx = xInverseBase->eval(b[i], xp);
    double xx = 0.0;
    auto ix = std::lower_bound(x.begin() + 1, x.end(), b[i]);
    if (ix == x.end())
    {
      xx = double(n);
    }
    else if (b[i] == *ix)
    {
      xx = double(std::distance(x.begin(), ix));
    }
    else
    {
      xx = double(std::distance(x.begin(), ix - 1)) + (b[i] - *(ix - 1)) / (*ix - *(ix - 1));
    }
    double xx1 = -cos(M_PI*xx / n);
    p[i] = base.eval(xx1, y);
  }
  return fitBase;
}

//=====================================================================================

  class TestChebfun : public Mantid::API::Algorithm {
  public:
    virtual const std::string name() const { return "TestChebfun"; }
    virtual int version() const { return 1; }
    virtual const std::string summary() const { return "Test new Chebfun functionality."; }
  private:
    void init();
    void exec();
  };

  /// Init this algorithm.
  void TestChebfun::init() {
    declareProperty(std::make_unique<Mantid::API::WorkspaceProperty<>>
      ("InputWorkspace", "", Mantid::Kernel::Direction::Input),
      "The input workspace");
    declareProperty("WorkspaceIndex", 0, "A workspace index");
    declareProperty("Error", 0.0, "An error value");
    declareProperty(
      std::make_unique<Mantid::API::WorkspaceProperty<>>("OutputWorkspace", "", Mantid::Kernel::Direction::Output),
      "Name of the output workspace");
  }

  /// Execute this algorithm.
  void TestChebfun::exec() {
    Mantid::API::MatrixWorkspace_sptr inputWorkspace = getProperty("InputWorkspace");
    auto nBins = inputWorkspace->blocksize();
    auto outputWorkspace = Mantid::API::WorkspaceFactory::Instance().create(inputWorkspace, 4, nBins, nBins);
    size_t wsIndex = static_cast<size_t>(int(getProperty("WorkspaceIndex")));


    setProperty("OutputWorkspace", outputWorkspace);
  }

  // Register the algorithm into the AlgorithmFactory
  DECLARE_ALGORITHM(TestChebfun)
}
}
}
