#ifndef MANTID_CURVEFITTING_COSTFUNCRWP_H_
#define MANTID_CURVEFITTING_COSTFUNCRWP_H_

//----------------------------------------------------------------------
// Includes
//----------------------------------------------------------------------
#include "MantidCurveFitting/CostFuncFitting.h"
#include "MantidCurveFitting/GSLMatrix.h"
#include "MantidCurveFitting/GSLVector.h"

namespace Mantid
{
namespace Kernel
{
  class Logger;
}
namespace CurveFitting
{
  class SeqDomain;
  class ParDomain;

/** Cost function for Rwp ... ... TODO: Make it complete!

    @author
    @date

    Copyright &copy;

    This file is part of Mantid.

    Mantid is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    Mantid is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    File change history is stored at: <https://github.com/mantidproject/mantid>.
    Code Documentation is available at: <http://doxygen.mantidproject.org>
*/
class DLLExport CostFuncRwp : public CostFuncFitting
{
public:
  /// Constructor
  CostFuncRwp();
  /// Virtual destructor
  virtual ~CostFuncRwp() {}

  /// Get name of minimizer
  virtual std::string name() const { return "Rwp";}

  /// Get short name of minimizer - useful for say labels in guis
  virtual std::string shortName() const {return "Rwp";}

  /// Calculate value of cost function
  virtual double val() const;

  /// Calculate the derivatives of the cost function
  /// @param der :: Container to output the derivatives
  virtual void deriv(std::vector<double>& der) const;

  /// Calculate the value and the derivatives of the cost function
  /// @param der :: Container to output the derivatives
  /// @return :: The value of the function
  virtual double valAndDeriv(std::vector<double>& der) const;

  virtual double valDerivHessian(bool evalFunction = true, bool evalDeriv = true, bool evalHessian = true) const;
  const GSLVector& getDeriv() const;
  const GSLMatrix& getHessian() const;
  void push();
  void pop();
  void drop();

  void setParameters(const GSLVector& params);
  void getParameters(GSLVector& params) const;

protected:

  virtual void calActiveCovarianceMatrix(GSLMatrix& covar, double epsrel = 1e-8);

  void addVal(
    API::FunctionDomain_sptr domain,
    API::FunctionValues_sptr values
    )const;
  void addValDerivHessian(
    API::IFunction_sptr function,
    API::FunctionDomain_sptr domain,
    API::FunctionValues_sptr values,
    bool evalFunction = true, bool evalDeriv = true, bool evalHessian = true) const;

private:

  mutable double m_value;
  mutable GSLVector m_der;
  mutable GSLMatrix m_hessian;

  mutable bool m_pushed;
  mutable double m_pushedValue;
  mutable GSLVector m_pushedParams;

  friend class SeqDomain;
  friend class ParDomain;

  Kernel::Logger & m_log;
};

} // namespace CurveFitting
} // namespace Mantid

#endif /*MANTID_CURVEFITTING_COSTFUNCRWP_H_*/
