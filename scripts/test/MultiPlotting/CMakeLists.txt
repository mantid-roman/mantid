#
## Tests for Muon GUIs
##

set ( TEST_PY_FILES
   AxisChangerTwoPresenter_test.py
   AxisChangerTwoView_test.py
   QuickEditPresenter_test.py
   QuickEditWidget_test.py
   MultiPlottingContext_test.py
   MultiPlotWidget_test.py
   Subplot_test.py
   SubPlotContext_test.py
)

check_tests_valid ( ${CMAKE_CURRENT_SOURCE_DIR} ${TEST_PY_FILES} )

# Prefix for test name=PythonAlgorithms
set ( PYUNITTEST_QT_API pyqt ) # force to use qt4
pyunittest_add_test ( ${CMAKE_CURRENT_SOURCE_DIR} python.MuonQt4 ${TEST_PY_FILES} )
unset ( PYUNITTEST_QT_API )
