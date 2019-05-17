// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//     NScD Oak Ridge National Laboratory, European Spallation Source
//     & Institut Laue - Langevin
// SPDX - License - Identifier: GPL - 3.0 +
// Mantid Coding standars <http://www.mantidproject.org/Coding_Standards>
#include "MantidCurveFitting/Functions/SqwFunction.h"
#include "MantidAPI/FunctionFactory.h"
#include "MantidAPI/MultiDomainFunction.h"
#include "MantidKernel/Logger.h"

namespace {
Mantid::Kernel::Logger g_log("SqwFunction");
}

namespace Mantid {
namespace CurveFitting {
namespace Functions {

DECLARE_FUNCTION(SqwFunction)

using namespace API;

/**
 * @brief Initialize elastic and inelastic parts, aliases, attributes, and ties
 */
//void DiffSphere::init() {
//  m_elastic = boost::dynamic_pointer_cast<ElasticDiffSphere>(
//      API::FunctionFactory::Instance().createFunction("ElasticDiffSphere"));
//  this->addFunction(m_elastic);
//  m_inelastic = boost::dynamic_pointer_cast<InelasticDiffSphere>(
//      API::FunctionFactory::Instance().createFunction("InelasticDiffSphere"));
//  this->addFunction(m_inelastic);
//
//  this->setAttributeValue("NumDeriv", true);
//  this->declareAttribute("Q", API::IFunction::Attribute(1.0));
//
//  // Set the aliases
//  this->setAlias("f1.Intensity", "Intensity");
//  this->setAlias("f1.Radius", "Radius");
//  this->setAlias("f1.Diffusion", "Diffusion");
//  this->setAlias("f1.Shift", "Shift");
//
//  // Set the ties between Elastic and Inelastic parameters
//  this->addDefaultTies(
//      "f0.Height=f1.Intensity,f0.Radius=f1.Radius,f0.Centre=0");
//  this->applyTies();
//}

SqwFunction::SqwFunction()
{
  declareAttribute("FWHMFormula", Attribute("d*x^2"));
  auto fwhmFunction = FunctionFactory::Instance().createFunction("UserFunction");
  fwhmFunction->setAttributeValue("Formula", "d*x^2");
  CompositeFunction::addFunction(fwhmFunction);
  CompositeFunction::addFunction(FunctionFactory::Instance().createFunction("MultiDomainFunction"));
}

size_t SqwFunction::addFunction(IFunction_sptr f) {
  if (m_useParent) return CompositeFunction::addFunction(f);
  auto mdFun = getInnerFunction();
  auto const i = mdFun->nFunctions();
  mdFun->addFunction(f);
  mdFun->setDomainIndex(i, i);
  m_useParent = true;
  checkFunction();
  m_useParent = false;
  return i;
}

IFunction_sptr SqwFunction::getFunction(std::size_t i) const {
  //if (m_useParent) return CompositeFunction::getFunction(i);
  return getInnerFunction()->getFunction(i);
}

std::size_t SqwFunction::nFunctions() const {
  //if (m_useParent) return CompositeFunction::nFunctions();
  return getInnerFunction()->nFunctions();
}

std::size_t SqwFunction::getNumberDomains() const {
  //if (m_useParent) return CompositeFunction::nFunctions();
  return getInnerFunction()->getNumberDomains();
}

std::vector<IFunction_sptr> SqwFunction::createEquivalentFunctions() const {
  return getInnerFunction()->createEquivalentFunctions();
}

//void SqwFunction::setParameter(const std::string &name, const double &value,
//  bool explicitlySet = true) {
//  
//}
//
//double SqwFunction::getParameter(const std::string &name) const {}
//
//bool SqwFunction::hasParameter(const std::string &name) const {}
//
//size_t SqwFunction::parameterIndex(const std::string &name) const {}

void SqwFunction::parseName(const std::string &varName, size_t &index, std::string &name) const {
}

MultiDomainFunction_sptr SqwFunction::getInnerFunction() const {
  m_useParent = true;
  auto fun = CompositeFunction::getFunction(1);
  m_useParent = false;
  return boost::dynamic_pointer_cast<MultiDomainFunction>(fun);
}

} // namespace Functions
} // namespace CurveFitting
} // namespace Mantid
