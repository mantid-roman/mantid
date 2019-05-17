// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//     NScD Oak Ridge National Laboratory, European Spallation Source
//     & Institut Laue - Langevin
// SPDX - License - Identifier: GPL - 3.0 +
#ifndef ENDERFCTEST_H_
#define ENDERFCTEST_H_

#include <cxxtest/TestSuite.h>

#include "MantidCurveFitting/Functions/SqwFunction.h"
#include "MantidAPI/FunctionFactory.h"
#include "MantidAPI/AlgorithmFactory.h"
#include "MantidAPI/Algorithm.h"

#include <iostream>

using namespace Mantid::API;
using namespace Mantid::CurveFitting::Functions;

class SqwFunctionTest : public CxxTest::TestSuite {
public:

  void test_values() {
    SqwFunction fun;
    fun.addFunction(FunctionFactory::Instance().createInitialized("name=Lorentzian"));
    fun.addFunction(FunctionFactory::Instance().createInitialized("name=Lorentzian"));
    fun.addFunction(FunctionFactory::Instance().createInitialized("name=Lorentzian"));
    
    fun.setParameter(3, 12);
    TS_ASSERT_DELTA(fun.getParameter(3), 12.0, 15);
    fun.setParameter(0, 1.2);
    TS_ASSERT_DELTA(fun.getParameter(0), 1.2, 15);
    //fun.setParameter("f0.d", 1.2);
    std::cerr << fun.getFunction(0)->asString() << std::endl << std::endl;
    std::cerr << fun.asString() << std::endl;
    TS_ASSERT_EQUALS(fun.nParams(), 10);
    TS_ASSERT_EQUALS(fun.nFunctions(), 3);
    TS_ASSERT_EQUALS(fun.getNumberDomains(), 3);



    auto loader = AlgorithmFactory::Instance().create("Load", -1);
    loader->initialize();
    loader->setProperty("Filename", "iris26176_graphite002_red.nxs");
    loader->setProperty("OutputWorkspace", "data");
    loader->execute();
  }
};

#endif /*ENDERFCTEST_H_*/
