#ifndef MANTID_CURVEFITTING_INTERVALCHEBFUN_H
#define MANTID_CURVEFITTING_INTERVALCHEBFUN_H

#include "MantidCurveFitting/Functions/SimpleChebfun.h"
#include <ostream>

namespace Mantid {
namespace CurveFitting {
namespace Functions {

typedef std::pair<double, double> Interval;

class IntervalChebfun : public Mantid::CurveFitting::Functions::SimpleChebfun
{
public:
  IntervalChebfun(const Mantid::CurveFitting::Functions::SimpleChebfun& fun);
  IntervalChebfun(const std::vector<double>& x, const std::vector<double>& y);
  Interval getInterval(const Interval& range) const;
  Interval getInterval(double x1, double x2) const { return getInterval(Interval(x1, x2)); }
private:
  void init();
  //Mantid::CurveFitting::SimpleChebfun m_deriv;
  std::vector<double> m_minima, m_minValues;
  std::vector<double> m_maxima, m_maxValues;
};

}
}
}

inline std::ostream& operator<<(std::ostream& ostr, const Mantid::CurveFitting::Functions::Interval& r)
{
  ostr << '(' << r.first << ',' << r.second << ')';
  return ostr;
}

#endif // MANTID_CURVEFITTING_INTERVALCHEBFUN_H