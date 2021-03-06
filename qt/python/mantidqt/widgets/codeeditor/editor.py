# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2017 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
#  This file is part of the mantidqt package
#
#
from __future__ import (absolute_import, unicode_literals)

# 3rd party imports

# local imports
from mantidqt.utils.qt import import_qt

# Import single-file editor from C++ wrapping
CodeEditor = import_qt('.._common', 'mantidqt.widgets', 'ScriptEditor')
