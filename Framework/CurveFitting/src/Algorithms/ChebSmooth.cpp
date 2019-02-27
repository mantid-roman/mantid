#include "MantidAPI/WorkspaceProperty.h"
#include "MantidAPI/WorkspaceFactory.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidCurveFitting/Functions/SimpleChebfun.h"
#include "MantidCurveFitting/Functions/TabulatedFunction.h"
#include "MantidCurveFitting/Functions/CubicSpline.h"
#include "MantidCurveFitting/Algorithms/Fit.h"

#include "MantidCurveFitting/Algorithms/ChebSmooth.h"
#include "MantidCurveFitting/Functions/IntervalChebfun.h"
#include "MantidCurveFitting/Functions/ChebfunUtils.h"

#include <boost/lexical_cast.hpp>
#include <numeric>
#include <iostream>

using namespace Mantid::CurveFitting;
using namespace Mantid::CurveFitting::Algorithms;
using namespace Mantid::CurveFitting::Functions;

// Register the algorithm into the AlgorithmFactory
DECLARE_ALGORITHM(ChebSmooth)

/// Init this algorithm.
void ChebSmooth::init() {
  declareProperty(
    std::make_unique<Mantid::API::WorkspaceProperty<>>("InputWorkspace", "", Mantid::Kernel::Direction::Input),
      "The input workspace");
  declareProperty("WorkspaceIndex", 0, "A workspace index");
  declareProperty("Error", 0.0, "An error value");
  declareProperty("DoAveraging", true, "Switch for averaging background points");
  declareProperty(
      std::make_unique<Mantid::API::WorkspaceProperty<>>("OutputWorkspace", "", Mantid::Kernel::Direction::Output),
      "Name of the output workspace");
  declareProperty(
      std::make_unique<Mantid::API::WorkspaceProperty<>>("AnalysisResult", "", Mantid::Kernel::Direction::Output),
      "Name of the test workspace");
}

/// Execute this algorithm.
void ChebSmooth::exec() {
  Mantid::API::MatrixWorkspace_sptr inputWorkspace = getProperty("InputWorkspace");
  bool doAveraging = getProperty("DoAveraging");

  auto nBins = inputWorkspace->blocksize();
  auto outputWorkspace = Mantid::API::WorkspaceFactory::Instance().create(inputWorkspace, 2, nBins, nBins);
  size_t wsIndex = static_cast<size_t>(int(getProperty("WorkspaceIndex")));
  setProperty("OutputWorkspace", outputWorkspace);

  m_progress = std::make_unique<Mantid::API::Progress>(this, 0.0, 1.0, 5);

  // Read the spectrum
  auto &X = inputWorkspace->readX(wsIndex);
  auto &Y = inputWorkspace->readY(wsIndex);

  //CHECK_OUT("x",X);
  //CHECK_OUT("y",Y);

  // Approximate the spectrum with a SimpleChebfun. It does some smoothing.
  SimpleChebfun cheb(X, Y);

  m_progress->report();
  std::cerr << "First smooth" << std::endl;

  // Evaluate the chebfun on the spectrum's x-points
  auto x = X;
  auto y = cheb(x);
  // Get the chebyshev expansion coefficients
  auto c = cheb.coeffs();

  // Save the smoothed spectrum and the cheb. coefficients to the first output spectrum
  outputWorkspace->dataX(0).assign(x.begin(), x.end());
  outputWorkspace->dataY(0).assign(y.begin(), y.end());
  outputWorkspace->dataE(0).assign(c.begin(), c.end());

  // Calculate the residual and save it in the second output spectrum
  std::transform(y.begin(), y.end(), Y.begin(), y.begin(), std::minus<double>());
  outputWorkspace->dataX(1).assign(x.begin(), x.end());
  outputWorkspace->dataY(1).assign(y.begin(), y.end());

  m_progress->report();

  // Estimate the noise
  std::transform(y.begin(), y.end(), y.begin(), [](double a){return a*a;});
  double sigma = getProperty("Error");
  if (sigma == 0.0)
  {
    sigma = std::accumulate(y.begin(), y.end(), 0.0) / nBins;
    sigma = sqrt(sigma) * 2;
  }
  std::cerr << "Sigma = " << sigma << std::endl;

  IntervalChebfun icheb(cheb);
  Interval range(cheb.startX(), cheb.endX());
  //std::cerr << icheb.getInterval(range) << std::endl;

  double sigma1 = 0.99 * sigma;
  double sigma2 = 1.05 * sigma;
  double dx = x[1] - x[0];
  const double xBegin = x.front();
  const double xEnd = x.back();
  double xp = xBegin;
  const double maxDx = (xEnd - xBegin) / 10;
  const double minDx = dx;
  std::vector<double> xx, yy, bg;
  xx.push_back(xp);
  yy.push_back(Y[0]);
  bg.push_back(Y[0]);

  double minLog = 1e100;
  double maxLog = 0.0;

  m_progress->report("Interval analysis");

  size_t count = 0;
  while(xp < xEnd) {
    double xp1 = xp + dx;
    if (xp1 >= xEnd) break;
    auto res = icheb.getInterval(xp, xp1);
    double diff = res.second - res.first;
    if (diff < sigma1 && dx < maxDx)
    {
      if (xp1 >= xEnd) break;
      //std::cerr << "diff " << xp << ' ' << xp1 << ' ' << diff << " increasing dx" << std::endl;
      dx *= 2;
    }
    else if (diff > sigma2 && dx > minDx)
    {
      //std::cerr << "diff " << xp << ' ' << xp1 << ' ' << diff << " decreasing dx" << std::endl;
      dx *= 0.75;
    }
    else
    {
      // y value at the end of current interval
      auto yp = cheb(xp1);
      // non-zero bp marks a point on the background
      auto bp = 0.0;

      auto n = xx.size();
      if (doAveraging && n > 2) {
        if (fabs(yp - yy[n-2]) < sigma1) {
          auto av = (yy[n-1] + yp) / 2;
          yy[n-1] = av;
          yp = av;
        }
      }

      auto dyb = fabs(yp - yy.back());
      if (dyb < sigma1) {
        //std::cerr << "Point outside dSigma " << xp1 << ' ' << dyb << ' ' << diff << std::endl;
        //bg.back() = yy.back();
        bp = yp;
      }

      xx.push_back(xp1);
      yy.push_back(yp);
      bg.push_back(bp);

      //std::cerr << xp << ' ' << xp1 << ' ' << yp << std::endl;

      xp = xp1;
      double ldx = log(dx);
      if (ldx < minLog) minLog = ldx;
      if (ldx > maxLog) maxLog = ldx;
    }
    interruption_point();
    ++count;
    //std::cerr << xp1 << ' ' << xEnd << ' ' << diff << ' ' << res << std::endl;
    //if (count > 1000) break;
  }
  if (xx.back() < xEnd)
  {
    xx.push_back(xEnd);
    yy.push_back(cheb(xEnd));
  }

  // ??!??!
  double midDx = exp(minLog + (maxLog - minLog)*0.4);

  std::vector<double> average(bg.size());
  size_t j = 1, nj = 0;
  double av = 0.0;
  double x1 = xx[j];
  double x0 = xx[j-1];
  double grad = (yy[j] - yy[j-1]) / (x1 - x0);
  for(size_t i = 0; i < Y.size(); ++i) {
    const double y = Y[i];
    if (X[i] < x1) {
      av += y - (yy[j-1] + grad * (X[i] - x0));
      ++nj;
    } else if (j < xx.size()) {
      ++j;
      x0 = x1;
      x1 = xx[j];
      grad = (yy[j] - yy[j-1]) / (x1 - x0);
      average[j] = av / (nj > 0 ? nj : 1);
      av = 0.0;
    } else {
      break;
    }
  }

  size_t iStart = 0;
  size_t iEnd = 0;
  const size_t iLast = bg.size() - 1;
  bool inBackground = false;
  for(size_t i = 2; i < bg.size(); ++i) {
    if (bg[i] != 0.0) {
      if (!inBackground) {
        iStart = i - 1;
        bg[iStart] = yy[iStart];
      }
      iEnd = i;
      inBackground = true;
    } else if (inBackground && iEnd < iLast) {
      inBackground = false;
      if (bg[iStart] > yy[iStart - 1] && bg[iEnd] > yy[iEnd + 1]) {
        std::cerr << "Bg: " << iStart << ' ' << iEnd << std::endl;
        std::fill(bg.begin() + iStart, bg.begin() + iEnd + 1, 0.0);
      } else {
        std::copy(yy.begin() + iStart, yy.begin() + iEnd + 1, bg.begin() + iStart);
      }
    }
  }

  auto testWorkspace = Mantid::API::WorkspaceFactory::Instance().create(inputWorkspace, 3, xx.size(), xx.size());
  setProperty("AnalysisResult", testWorkspace);

  testWorkspace->dataX(0).assign(xx.begin(), xx.end());
  testWorkspace->dataY(0).assign(yy.begin(), yy.end());

  testWorkspace->dataX(1).assign(xx.begin(), xx.end());
  testWorkspace->dataY(1).assign(bg.begin(), bg.end());

  testWorkspace->dataX(2).assign(xx.begin(), xx.end());
  testWorkspace->dataY(2).assign(average.begin(), average.end());

  std::cerr << "Smaller size " << xx.size() << std::endl;

  //CHECK_OUT("xx",xx);
  //CHECK_OUT("yy",yy);

  m_progress->report();

}
