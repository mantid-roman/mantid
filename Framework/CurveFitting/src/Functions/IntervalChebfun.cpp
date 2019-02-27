#include "MantidCurveFitting/Functions//IntervalChebfun.h"

#include <iostream>

namespace Mantid {
namespace CurveFitting {
namespace Functions {

IntervalChebfun::IntervalChebfun(const SimpleChebfun& fun) :
  SimpleChebfun(fun.order(), fun, fun.startX(), fun.endX())
{
  init();
}

IntervalChebfun::IntervalChebfun(const std::vector<double>& x, const std::vector<double>& y) :
  SimpleChebfun(x, y)
{
  init();
}

void IntervalChebfun::init()
{
  auto der1 = derivative();
  auto roots = der1.roughRoots();
  auto der2 = der1.derivative();
  for (auto r = roots.begin(); r != roots.end(); ++r)
  {
    auto v = (*this)(*r);
    if (der2(*r) > 0.0)
    {
      m_minima.push_back(*r);
      m_minValues.push_back(v);
      //std::cerr << "Minimum at " << *r << " is " << v << std::endl;
    }
    else
    {
      m_maxima.push_back(*r);
      m_maxValues.push_back(v);
      //std::cerr << "Maximum at " << *r << " is " << v << std::endl;
    }
  }
}

Interval IntervalChebfun::getInterval(const Interval& range) const
{
  double v1 = (*this)(range.first);
  double v2 = (*this)(range.second);
  double minVal = std::min(v1, v2);
  auto minBegin = std::lower_bound(m_minima.begin(), m_minima.end(), range.first);
  auto minEnd = std::upper_bound(minBegin, m_minima.end(), range.second);
  for (auto it = minBegin; it != minEnd; ++it)
  {
    auto i = std::distance(m_minima.begin(), it);
    auto v = m_minValues[i];
    if (v < minVal) minVal = v;
  }
  double maxVal = std::max(v1, v2);
  auto maxBegin = std::lower_bound(m_maxima.begin(), m_maxima.end(), range.first);
  auto maxEnd = std::upper_bound(maxBegin, m_maxima.end(), range.second);
  for (auto it = maxBegin; it != maxEnd; ++it)
  {
    auto i = std::distance(m_maxima.begin(), it);
    auto v = m_maxValues[i];
    if (v > maxVal) maxVal = v;
  }
  return Interval(minVal, maxVal);
}

}
}
}