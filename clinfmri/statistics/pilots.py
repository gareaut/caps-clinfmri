#! /usr/bin/env python
##########################################################################
# CAPS - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import sys
import unittest
import os
from PySide import QtGui
import logging
import datetime

logging.basicConfig(level=logging.INFO)

# CAPSUL import
from capsul.qt_gui.widgets import PipelineDevelopperView
from capsul.study_config.study_config import StudyConfig
from capsul.api import get_process_instance


pipeline = get_process_instance("clinfmri.statistics.spm_first_level_pipeline")
pipeline.smoother_switch = "no_smoothing"
pipeline.complete_regressors = "yes"
app = QtGui.QApplication(sys.argv)
view1 = PipelineDevelopperView(pipeline, allow_open_controller=True, show_sub_pipelines=True)
view1.show()
app.exec_()
del view1

