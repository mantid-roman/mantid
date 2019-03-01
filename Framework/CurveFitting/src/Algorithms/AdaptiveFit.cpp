// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//     NScD Oak Ridge National Laboratory, European Spallation Source
//     & Institut Laue - Langevin
// SPDX - License - Identifier: GPL - 3.0 +
#include "MantidCurveFitting/Algorithms/AdaptiveFit.h"
#include "MantidCurveFitting/Algorithms/Fit.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidAPI/Workspace.h"
#include <iostream>


namespace Mantid {
namespace CurveFitting {
namespace Algorithms {

using namespace Mantid::Kernel;
using namespace Mantid::API;

// Register the algorithm into the AlgorithmFactory
DECLARE_ALGORITHM(AdaptiveFit)

//----------------------------------------------------------------------------------------------

/// Algorithms name for identification. @see Algorithm::name
const std::string AdaptiveFit::name() const { return "AdaptiveFit"; }

/// Algorithm's version for identification. @see Algorithm::version
int AdaptiveFit::version() const { return 1; }

/// Algorithm's summary for use in the GUI and help. @see Algorithm::summary
const std::string AdaptiveFit::summary() const {
  return "Fit a function on a workspace.";
}

//----------------------------------------------------------------------------------------------
/// Initialize the algorithm's properties.
void AdaptiveFit::initConcrete() {
  declareProperty(make_unique<WorkspaceProperty<API::MatrixWorkspace>>(
                      "OutputWorkspace", "", Direction::Output),
                  "An output workspace.");
}

//----------------------------------------------------------------------------------------------
/// Execute the algorithm.
void AdaptiveFit::execConcrete() {

  Workspace_sptr inputWS = getProperty("InputWorkspace");
  Fit fit;
  fit.setChild(true);
  fit.initialize();
  fit.setProperty("Function", m_function);
  fit.setProperty("InputWorkspace", inputWS);
  fit.setProperty("CreateOutput", true);
  fit.execute();
  bool ok = fit.isExecuted();
  MatrixWorkspace_sptr outputWS = fit.getProperty("OutputWorkspace");
  
  // Store the result.
  setProperty("OutputWorkspace", outputWS);
}

} // namespace Algorithms
} // namespace CurveFitting
} // namespace Mantid
