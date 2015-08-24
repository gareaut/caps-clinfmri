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
from capsul.process.loader import get_process_instance

# CAPS import
from caps.toy_datasets import get_sample_data


# Configure the environment
start_time = datetime.datetime.now()
print "Start Configuration", start_time
study_config = StudyConfig(
    modules=["MatlabConfig", "SPMConfig", "NipypeConfig", "FSLConfig",
             "FreeSurferConfig", "SmartCachingConfig"],
    matlab_exec="/neurospin/local/bin/matlab",
    spm_directory="/i2bm/local/spm8-6313",
    use_matlab=True,
    use_spm=True,
    use_nipype=True,
    use_smart_caching=True,
    output_directory="/volatile/nsap/catalogue/spm_first_level/")
print "Done in {0} seconds".format(datetime.datetime.now() - start_time)


# Create pipeline
start_time = datetime.datetime.now()
print "Start Pipeline Creation", start_time
pipeline = get_process_instance(
    "caps.nsap.functional_statistic.pipeline.spm_first_level_pipeline.xml")
print "Done in {0} seconds.".format(datetime.datetime.now() - start_time)


# Set pipeline input parameters
start_time = datetime.datetime.now()
print "Start Parametrization", start_time
localizer_dataset = get_sample_data("localizer")
pipeline.behavioral_data = [localizer_dataset.onsets]
pipeline.fmri_sessions = [localizer_dataset.preproc_fmri]
pipeline.time_repetition = localizer_dataset.TR
pipeline.realignment_parameters = localizer_dataset.mouvment_parameters
pipeline.condition_name = "Conditions"
pipeline.onset_name = "Onsets"
pipeline.duration_name = "Durations"
pipeline.delimiter = ";"
pipeline.start = 0
pipeline.contrasts = [
    ("Horizontal Checkerboard","T",['damier_H',],[1,]),
    ("Vertical Checkerboard","T",['damier_V',],[1,]),
    ("Right Audio Click","T",['clicDaudio',],[1,]),
    ("Left Audio Click","T",['clicGaudio',],[1,]),
    ("Right Video Click","T",['clicDvideo',],[1,]),
    ("Left Video Click","T",['clicGvideo',],[1,]),
    ("Audio Computation","T",['calculaudio',],[1,]),
    ("Video Computation","T",['calculvideo',],[1,]),
    ("Video Sentences","T",['phrasevideo',],[1,]),
    ("Audio Sentences","T",['phraseaudio',],[1,]),
    ("Left Click","T",['clicGaudio','clicGvideo'],[0.5,0.5]),
    ("Right Click","T",['clicDaudio','clicDvideo'],[0.5,0.5]),
    ("Checkerboard","T",['damier_H','damier_V'],[0.5,0.5]),
    ("Motor","T",['clicGaudio','clicGvideo','clicDaudio','clicDvideo'],[0.25,0.25,0.25,0.25]),
    ("Computation","T",['calculaudio','calculvideo'],[0.5,0.5]),
    ("Sentences","T",['phrasevideo','phraseaudio'],[0.5,0.5]),
    ("Audio","T",['clicDaudio','clicGaudio','calculaudio','phraseaudio'],[0.25,0.25,0.25,0.25]),
    ("Video","T",['clicDvideo','clicGvideo','calculvideo','phrasevideo'],[0.25,0.25,0.25,0.25]),
    ("H Checkerboard - V Checkerboard","T",['damier_H','damier_V'],[1,-1]),
    ("V Checkerboard - H Checkerboard","T",['damier_H','damier_V'],[-1,1]),
    ("Left Click - Right Click","T",
        ['clicGaudio','clicGvideo','clicDaudio','clicDvideo'],[0.5,0.5,-0.5,-0.5]),
    ("Right Click - Left Click","T",
        ['clicGaudio','clicGvideo','clicDaudio','clicDvideo'],[-0.5,-0.5,0.5,0.5]),
    ("Video - Audio","T",
        ['clicDvideo','clicGvideo','calculvideo','phrasevideo','clicDaudio','clicGaudio',
        'calculaudio','phraseaudio'],[0.25,0.25,0.25,0.25,-0.25,-0.25,-0.25,-0.25]),
    ("Audio - Video","T",
        ['clicDvideo','clicGvideo','calculvideo','phrasevideo','clicDaudio','clicGaudio',
        'calculaudio','phraseaudio'],[-0.25,-0.25,-0.25,-0.25,0.25,0.25,0.25,0.25]),
    ("Motor - Cognitive","T",
        ['clicDaudio','clicGaudio','clicDvideo','clicGvideo','calculaudio','calculvideo',
        'phrasevideo','phraseaudio'],[0.25,0.25,0.25,0.25,-0.25,-0.25,-0.25,-0.25]),
    ("Cognitive - Motor","T",
        ['clicDaudio','clicGaudio','clicDvideo','clicGvideo','calculaudio','calculvideo',
        'phrasevideo','phraseaudio'],[-0.25,-0.25,-0.25,-0.25,0.25,0.25,0.25,0.25]),
    ("Audio Computation - Audio Sentences","T",['calculaudio','phraseaudio'],[1,-1]),
    ("Video Computation - Video Sentences","T",['calculvideo','phrasevideo'],[1,-1]),
    ("Computation - Sentences","T",
        ['calculaudio','calculvideo','phraseaudio','phrasevideo'],[0.5,0.5,-0.5,-0.5]),
    ("Video - Checkerboard","T",
        ['clicDvideo','clicGvideo','calculvideo','phrasevideo','damier_H','damier_V'],
        [0.25,0.25,0.25,0.25,-0.5,-0.5]),
    ("Video Sentences - Checkerboard","T",
        ['phrasevideo','damier_H','damier_V'],[1,-0.5,-0.5]),
    ("Audio Click - Audio Sentences","T",
        ['clicDaudio','clicGaudio','phraseaudio'],[0.5,0.5,-1]),
    ("Video Click - Video Sentences","T",
        ['clicDvideo','clicGvideo','phrasevideo'],[0.5,0.5,-1]),
]
print "Done in {0} seconds.".format(datetime.datetime.now() - start_time)


# View pipeline
if 0:
    app = QtGui.QApplication(sys.argv)
    view1 = PipelineDevelopperView(pipeline)
    view1.show()
    app.exec_()
    del view1

# Execute the pipeline in the configured study
study_config.run(pipeline, verbose=1)
