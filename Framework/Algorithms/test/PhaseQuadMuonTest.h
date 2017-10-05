#ifndef MANTID_ALGORITHMS_PHASEQUADMUONTEST_H_
#define MANTID_ALGORITHMS_PHASEQUADMUONTEST_H_

#include <math.h>

#include <cxxtest/TestSuite.h>
#include "MantidAPI/AlgorithmManager.h"
#include "MantidAPI/AnalysisDataService.h"
#include "MantidAPI/Axis.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidDataObjects/TableWorkspace.h"
#include "MantidAPI/TableRow.h"



using namespace Mantid::DataObjects;
using namespace Mantid::API;

namespace {

void populatePhaseTable(ITableWorkspace_sptr phaseTable,
                        std::vector<std::string> names) {
  phaseTable->addColumn("int", names[0]);
  phaseTable->addColumn("double", names[1]);
  phaseTable->addColumn("double", names[2]);

  for (int i = 0; i < 16; i++) {
    TableRow phaseRow1 = phaseTable->appendRow();
    phaseRow1 << i << 1. << 2.;
    TableRow phaseRow2 = phaseTable->appendRow();
    phaseRow2 << i << 1. << 2.;
  }
}
void populatePhaseTable(ITableWorkspace_sptr phaseTable) {
  populatePhaseTable(phaseTable, {"DetectorID", "Asymmetry", "Phase"});
}

IAlgorithm_sptr setupAlg(MatrixWorkspace_sptr m_loadedData, bool isChildAlg,
                         ITableWorkspace_sptr phaseTable) {
  // Set up PhaseQuad
  IAlgorithm_sptr phaseQuad = AlgorithmManager::Instance().create("PhaseQuad");
  phaseQuad->setChild(isChildAlg);
  phaseQuad->initialize();
  phaseQuad->setProperty("InputWorkspace", m_loadedData);
  phaseQuad->setProperty("PhaseTable", phaseTable);
  phaseQuad->setPropertyValue("OutputWorkspace", "outputWs");
  return phaseQuad;
}

IAlgorithm_sptr setupAlg(MatrixWorkspace_sptr m_loadedData, bool isChildAlg) {
  // Create and populate a detector table
  boost::shared_ptr<ITableWorkspace> phaseTable(
      new Mantid::DataObjects::TableWorkspace);
  populatePhaseTable(phaseTable);
  auto stuff =phaseTable->getColumnNames();
  return setupAlg(m_loadedData, isChildAlg, phaseTable);
}

IAlgorithm_sptr setupAlg(MatrixWorkspace_sptr m_loadedData, bool isChildAlg,
                         std::vector<std::string> names) {
  // Create and populate a detector table
  boost::shared_ptr<ITableWorkspace> phaseTable(
      new Mantid::DataObjects::TableWorkspace);
  populatePhaseTable(phaseTable, names);

  return setupAlg(m_loadedData, isChildAlg, phaseTable);
}

MatrixWorkspace_sptr loadMuonDataset() {
  IAlgorithm_sptr loader = AlgorithmManager::Instance().create("Load");
  loader->setChild(true);
  loader->initialize();
  loader->setProperty("Filename", "emu00006473.nxs");
  loader->setPropertyValue("OutputWorkspace", "outputWs");
  loader->execute();
  Workspace_sptr temp = loader->getProperty("OutputWorkspace");
  MatrixWorkspace_sptr m_loadedData =
      boost::dynamic_pointer_cast<MatrixWorkspace>(temp);
  return m_loadedData;
}
}

class PhaseQuadMuonTest : public CxxTest::TestSuite {
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static PhaseQuadMuonTest *createSuite() { return new PhaseQuadMuonTest(); }

  static void destroySuite(PhaseQuadMuonTest *suite) { delete suite; }

  void setUp() {
    if (!m_loadedData) {
      m_loadedData = loadMuonDataset();
    }
  }

  void testTheBasics() {
    IAlgorithm_sptr phaseQuad =
        AlgorithmManager::Instance().create("PhaseQuad");
    TS_ASSERT_EQUALS(phaseQuad->name(), "PhaseQuad");
    TS_ASSERT_EQUALS(phaseQuad->category(), "Muon");
    TS_ASSERT_THROWS_NOTHING(phaseQuad->initialize());
    TS_ASSERT(phaseQuad->isInitialized());

  }

  void testExecPhaseTable() {
    IAlgorithm_sptr phaseQuad = setupAlg(m_loadedData, true);
    TS_ASSERT_THROWS_NOTHING(phaseQuad->execute());
    TS_ASSERT(phaseQuad->isExecuted());

    // Get the output ws
    MatrixWorkspace_sptr outputWs = phaseQuad->getProperty("OutputWorkspace");

    TS_ASSERT_EQUALS(outputWs->getNumberHistograms(), 2);
    TS_ASSERT_EQUALS(
        outputWs->getSpectrum(0).readX(),
        m_loadedData->getSpectrum(0).readX()); // Check outputWs X values
    TS_ASSERT_EQUALS(outputWs->getSpectrum(1).readX(),
                     m_loadedData->getSpectrum(1).readX());

    const auto specReY = outputWs->getSpectrum(0).y();
    const auto specReE = outputWs->getSpectrum(0).e();
    const auto specImY = outputWs->getSpectrum(1).y();
    const auto specImE = outputWs->getSpectrum(1).e();
    // Check real Y values
    TS_ASSERT_DELTA(specReY[0], 2.1969, 0.0001);
    TS_ASSERT_DELTA(specReY[20], 0.0510, 0.0001);
    TS_ASSERT_DELTA(specReY[50], -0.0525, 0.0001);
    // Check real E values
    TS_ASSERT_DELTA(specReE[0], 0.0024, 0.0001);
    TS_ASSERT_DELTA(specReE[20], 0.0041, 0.0001);
    TS_ASSERT_DELTA(specReE[50], 0.0047, 0.0001);
    // Check imaginary Y values
    TS_ASSERT_DELTA(specImY[0], -0.1035, 0.0001);
    TS_ASSERT_DELTA(specImY[20], -0.0006, 0.0001);
    TS_ASSERT_DELTA(specImY[50], 0.0047, 0.0001);
    // Check imaginary E values
    TS_ASSERT_DELTA(specImE[0], 0.0002, 0.0001);
    TS_ASSERT_DELTA(specImE[20], 0.0004, 0.0001);
    TS_ASSERT_DELTA(specImE[50], 0.0005, 0.0001);
  }
  void testNoPhase() {
    std::vector<std::string> names = {"ID", "Asym", "dummy"};
    IAlgorithm_sptr phaseQuad = setupAlg(m_loadedData, true, names);
    TS_ASSERT_THROWS(phaseQuad->execute(), std::runtime_error);
  }
  void testNoAsymm() {
    std::vector<std::string> names = {"ID", "AsYMg", "phase"};
    MatrixWorkspace_sptr m_loadedData = loadMuonDataset();
    IAlgorithm_sptr phaseQuad = setupAlg(m_loadedData, true, names);
    TS_ASSERT_THROWS(phaseQuad->execute(), std::runtime_error);
  }
  void testTwoPhases() {
    std::vector<std::string> names = {"ID", "Phase", "phi"};
    IAlgorithm_sptr phaseQuad = setupAlg(m_loadedData, true, names);
    TS_ASSERT_THROWS(phaseQuad->execute(), std::runtime_error);
  }
  void testTwoAsymm() {
    std::vector<std::string> names = {"ID", "Asym", "Asymm"};
    IAlgorithm_sptr phaseQuad = setupAlg(m_loadedData, true, names);
    TS_ASSERT_THROWS(phaseQuad->execute(), std::runtime_error);
  }
  void testSwapOrder() {
    std::vector<std::string> names = {"ID", "phase", "Asymm"};
    IAlgorithm_sptr phaseQuad = setupAlg(m_loadedData, true, names);
    TS_ASSERT_THROWS_NOTHING(phaseQuad->execute());
    TS_ASSERT(phaseQuad->isExecuted());

    // Get the output ws
    MatrixWorkspace_sptr outputWs = phaseQuad->getProperty("OutputWorkspace");

    TS_ASSERT_EQUALS(outputWs->getNumberHistograms(), 2);
    TS_ASSERT_EQUALS(
        outputWs->getSpectrum(0).readX(),
        m_loadedData->getSpectrum(0).readX()); // Check outputWs X values
    TS_ASSERT_EQUALS(outputWs->getSpectrum(1).readX(),
                     m_loadedData->getSpectrum(1).readX());

    const auto specReY = outputWs->getSpectrum(0).y();
    const auto specReE = outputWs->getSpectrum(0).e();
    const auto specImY = outputWs->getSpectrum(1).y();
    const auto specImE = outputWs->getSpectrum(1).e();
    // Check real Y values
    TS_ASSERT_DELTA(specReY[0], -0.9358, 0.0001);
    TS_ASSERT_DELTA(specReY[20], -0.0236, 0.0001);
    TS_ASSERT_DELTA(specReY[50], 0.0197, 0.0001);
    // Check real E values
    TS_ASSERT_DELTA(specReE[0], 0.0009, 0.0001);
    TS_ASSERT_DELTA(specReE[20], 0.0018, 0.0001);
    TS_ASSERT_DELTA(specReE[50], 0.0020, 0.0001);
    // Check imaginary Y values
    TS_ASSERT_DELTA(specImY[0], -0.6379, 0.0001);
    TS_ASSERT_DELTA(specImY[20], -0.0129, 0.0001);
    TS_ASSERT_DELTA(specImY[50], 0.0159, 0.0001);
    // Check imaginary E values
    TS_ASSERT_DELTA(specImE[0], 0.0006, 0.0001);
    TS_ASSERT_DELTA(specImE[20], 0.0012, 0.0001);
    TS_ASSERT_DELTA(specImE[50], 0.0014, 0.0001);
  }
  // add test for different order

private:
  MatrixWorkspace_sptr m_loadedData;
};

class PhaseQuadMuonTestPerformance : public CxxTest::TestSuite {

public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static PhaseQuadMuonTestPerformance *createSuite() {
    return new PhaseQuadMuonTestPerformance();
  }

  static void destroySuite(PhaseQuadMuonTestPerformance *suite) {
    delete suite;
  }

  void setUp() override {
    m_loadedData = loadMuonDataset();
    phaseQuad = setupAlg(m_loadedData, false);
  }

  void tearDown() override {
    Mantid::API::AnalysisDataService::Instance().remove("outputWs");
  }

  void testPerformanceWs() { phaseQuad->execute(); }

private:
  MatrixWorkspace_sptr m_loadedData;
  IAlgorithm_sptr phaseQuad;
};

#endif /* MANTID_ALGORITHMS_PHASEQUADMUONTEST_H_ */
