#! /usr/bin/env python
##########################################################################
# CAPS - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import os
from PySide import QtGui
import datetime
from glob import glob
import shutil
import sys

#for item in sys.path:
    #if any(word in item for word in ['mmutils', 'clinfmri']):
        #sys.path.remove(item)
#sys.path.append('/volatile/git/caps-clinfmri')
#sys.path.append('/volatile/git/caps-mmutils')
import clinfmri
import mmutils
print clinfmri.__file__
print mmutils.__file__

# CAPSUL import
from capsul.qt_gui.widgets import PipelineDevelopperView
from capsul.study_config.study_config import StudyConfig

from capsul.process.loader import get_process_instance

# Configure the environment
start_time = datetime.datetime.now()

# create output_dir
out_dir = "/neurospin/brainomics/2015_TEMPO_imagen/capsul_preproc_out"
if os.path.isdir(out_dir):
    shutil.rmtree(out_dir)
os.makedirs(out_dir)

print "Start Configuration", start_time
study_config = StudyConfig(
    modules=["MatlabConfig", "FSLConfig", "SPMConfig", "NipypeConfig"],
    use_matlab=True,
    matlab_exec="/neurospin/local/bin/matlab",
    # spm_directory="/i2bm/local/spm8-6313",
    spm_directory="/neurospin/imagen/workspace/spm_libs/spm8-6313",
    fsl_config="/etc/fsl/4.1/fsl.sh",
    use_spm=True,
    use_nipype=True,
    use_fsl=True,
    number_of_cpus=1,
    generate_logging=True,
    use_scheduler=True,
    output_directory=out_dir)

print "Done in {0} seconds".format(datetime.datetime.now() - start_time)
print ""

# inputs parameters
fmri_sessions = []
behavioral_data = []
contrasts = []
realignment_parameters = []

nii_face_list = glob("/neurospin/brainomics/2015_TEMPO_imagen/"
                     "raw_BL_niftis/*/Session*/EPI*faces*/*.nii.gz")
nii_t1_list = glob("/neurospin/brainomics/2015_TEMPO_imagen/"
                   "raw_BL_niftis/*/Session*/ADNI_MPRAGE/00*.nii.gz")

# XXX only one subject
nii_face_list = nii_face_list[:1]
nii_t1_list = nii_t1_list[:1]

print nii_face_list
print nii_t1_list

# Create pipeline
start_time = datetime.datetime.now()
print "Start Pipeline Creation", start_time
pipeline = get_process_instance("clinfmri.preproc.fmri_preproc.xml")
print "Done in {0} seconds.".format(datetime.datetime.now() - start_time)

# View pipeline
if 0:
    pipeline_visu = get_process_instance("clinfmri.preproc.fmri_preproc_bbox.xml")
    app = QtGui.QApplication(sys.argv)
    view1 = PipelineDevelopperView(pipeline_visu)
    view1.show()
    app.exec_()

# now parameterize pipeline
print "Start Pipeline Parametrization", start_time
pipeline.fmri_file = nii_face_list
pipeline.structural_file = nii_t1_list
# configure slice timing
pipeline.select_slicer = "spm"
pipeline.force_repetition_time = 2.2
pipeline.force_slice_order = [40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29,
                              28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17,
                              16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4,
                              3, 2, 1]
pipeline.reorient_and_register = "do_nothing"
pipeline.ref_slice = 40
# configure realignment
pipeline.wrap = [0, 1, 0]
pipeline.write_wrap = [0, 1, 0]
pipeline.write_which = [1, 1]
pipeline.realign_register_to_mean = False
# config structural segment
pipeline.input_image_info = (0.0001, 60, (True, True))
pipeline.tissues = [(("/i2bm/local/spm8-6313/toolbox/Seg/TPM.nii", 1), 2,
                     (True, True),
                     (False, True)),
                    (("/i2bm/local/spm8-6313/toolbox/Seg/TPM.nii", 2), 2,
                     (True, True),
                     (False, True)),
                    (("/i2bm/local/spm8-6313/toolbox/Seg/TPM.nii", 3), 2,
                     (True, False),
                     (False, False)),
                    (("/i2bm/local/spm8-6313/toolbox/Seg/TPM.nii", 4), 3,
                     (False, False),
                     (False, False)),
                    (("/i2bm/local/spm8-6313/toolbox/Seg/TPM.nii", 5), 4,
                     (False, False),
                     (False, False)),
                    (("/i2bm/local/spm8-6313/toolbox/Seg/TPM.nii", 6), 2,
                     (False, False),
                     (False, False))
                    ]
pipeline.warping_regularization = 4
pipeline.affine_regularization = "mni"
pipeline.sampling_distance = 3
pipeline.write_deformation_fields = [False, True]
pipeline.mni_struct_template = ("/neurospin/imagen/src/scripts/FU2/processed"
                                "/MNI152_T1_1mm_brain_mask.nii")
# configure normalize
pipeline.select_registration = "template"
# pipeline.template_type = "3D"
pipeline.select_normalization = "fmri"
# pipeline.template_file = ("/neurospin/imagen/src/scripts/FU2/processed/"
#                          "mni_icbm152_t1_tal_nlin_asym_09c.nii")
# pipeline.template_type = "4D"
pipeline.template_file = ('/neurospin/imagen/workspace/fmri/scripts/'
			  'ImagenEPI200_3mm.nii')
print "Done in {0} seconds.".format(
    datetime.datetime.now() - start_time)

# run pipeline
study_config.run(pipeline, verbose=1)
