#ifndef MANTID_CURVEFITTING_CHEBSMOOTH_H_
#define MANTID_CURVEFITTING_CHEBSMOOTH_H_

#include "MantidAPI/Algorithm.h"

namespace Mantid {
namespace CurveFitting {
namespace Algorithms {

class ChebSmooth : public Mantid::API::Algorithm {
public:
  virtual const std::string name() const { return "ChebSmooth"; }
  virtual int version() const { return 1; }
  virtual const std::string summary() const { return "Smooth a spectrum using Chebfun."; }
private:
  void init();
  void exec();
  std::unique_ptr<Mantid::API::Progress> m_progress = nullptr;
};

}
}
}

#endif // MANTID_CURVEFITTING_CHEBSMOOTH_H_
