// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2007 ISIS Rutherford Appleton Laboratory UKRI,
//     NScD Oak Ridge National Laboratory, European Spallation Source
//     & Institut Laue - Langevin
// SPDX - License - Identifier: GPL - 3.0 +
#ifndef MANTID_SQWFUNCTION_H_
#define MANTID_SQWFUNCTION_H_

#include "MantidAPI/IFunction_fwd.h"
#include "MantidAPI/ImmutableCompositeFunction.h"

namespace Mantid {
namespace CurveFitting {
namespace Functions {

using namespace API;

class DLLExport SqwFunction : public API::ImmutableCompositeFunction {

public:
  SqwFunction();
  std::string name() const override { return "SqwFunction"; }
  const std::string category() const override { return "QuasiElastic"; }
  size_t addFunction(IFunction_sptr f) override;
  IFunction_sptr getFunction(std::size_t i) const override;
  std::size_t nFunctions() const override;
  std::size_t getNumberDomains() const override;
  std::vector<IFunction_sptr> createEquivalentFunctions() const override;

  //void setParameter(const std::string &name, const double &value,
  //  bool explicitlySet = true) override;
  ////void setParameterDescription(const std::string &name,
  ////  const std::string &description) override;
  //double getParameter(const std::string &name) const override;
  //bool hasParameter(const std::string &name) const override;
  //size_t parameterIndex(const std::string &name) const override;

private:
  void parseName(const std::string &varName, size_t &index, std::string &name) const override;
  MultiDomainFunction_sptr getInnerFunction() const;
  mutable bool m_useParent = false;
};

} // namespace Functions
} // namespace CurveFitting
} // namespace Mantid

#endif /*MANTID_SQWFUNCTION_H_*/
