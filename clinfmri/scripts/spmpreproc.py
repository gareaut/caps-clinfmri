#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013-2016
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
from __future__ import print_function
import argparse
import os
import shutil
import numpy

# Bredala import
try:
    import bredala
    bredala.USE_PROFILER = False
    bredala.register("qap.viz.plotting", names=["plot_mosaic"])
    bredala.register("capsul.process", names=["get_process_instance"])
    bredala.register("mmutils.adapters.io", names=["gzip_file"])
except:
    pass

# Mmutils import
from mmutils.adapters.io import gzip_file

# CAPSUL import
from capsul.study_config import StudyConfig
from capsul.process import get_process_instance

# Qap import
from qap.viz.plotting import plot_mosaic


# Parameters to keep trace
__hopla__ = ["sid", "tmp_capsul", "outdir", "funcfile", "logdir"]


# Script documentation
doc = """
fMRI spatial preprocessings
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Based on SPM and FSL realize fMRI sequences preprocessings. It is
possible to select the slice timing and the normalization algorithms
by setting the 'select_slicer' and 'select_normalization' pipeline 
parameters. We can choose to use the 'fsl' or 'spm' slicer or we
can select a 'fmri' or 'anat' normalization:

* in the case of an 'anat' normalization a t1 image has to be
  specified in the 'structural_file' parameter. This case simply
  register the t1 image to the mean fMRI image and then register the t1
  image in the fMRI space to the MNI template using the unified
  'New Segment' SPM procedure.
* in the cas of a 'fmri' normalization a functional template image has
  to be specified in the 'template_file'. This case simply register
  the mean fMRI image to the template.

In both case a 'fmri_file' containing the functional sequence has to be
specified. The latter is expected to contain the repetition time and
slice orders in this header. If it is not the case the
'force_slice_orders' and 'force_repetition_time' parameters have to be
specified.

The generated data are organized with respect to the openfmri standard.

Softawares involved:

* SPM
* FSL: for the slice timing as it is designed for multi-band
  acquisitions and the brain extraction.

Steps (with anatomical template alignement):
    * Slice timing: correct differences in image acquisition time
      between slices.
    * Realignement: estimate within modality volume to volume rigid
      body alignment.
    * Coregistration: register the t1 image to the functional space.
    * Normalization: une the SPM unified 'New Segment' algorithm to
      register the t1 image in the functional space to the MNI
      template space.
    * Smoothing
    * BET: brain extraction in the functional sequence.

Steps (with functional template alignement):
    * Slice timing: correct differences in image acquisition time
      between slices.
    * Realignement: estimate within modality volume to volume rigid
      body alignment.
    * Normalization: register the mean functional image to the template
      space.
    * Smoothing
    * BET: brain extraction in the functional sequence.

Command:

python $VIRTUAL_ENV/git/caps-clinfmri/clinfmri/scripts/spmpreproc.py \
    -v 2 \
    -c 000000001274 \
    -f /etc/fsl/4.1/fsl.sh \
    -s /i2bm/local/bin/spm12 \
    -o /volatile/nsap/spmpreproc \
    -n fmri \
    -t spm \
    -u 2.2 \
    -l 40 39 38 37 36 35 34 33 32 31 30 29 28 27 26 25 24 23 22 21 20 19 18 \
       17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 \
    -i /neurospin/imagen/BL/processed/nifti/000000001274/SessionA/EPI_stop_signal/000000001274s701a1007.nii.gz \
    -r /neurospin/imagen/workspace/fmri/scripts/ImagenEPI200_3mm.nii \
    -k /neurospin/imagen/workspace/fmri/scripts/ImagenEPI200_3mm_brain_mask.nii.gz \
    -w 40 \
    -e 
"""


def is_file(filearg):
    """ Type for argparse - checks that file exists but does not open.
    """
    if not os.path.isfile(filearg):
        raise argparse.ArgumentError(
            "The file '{0}' does not exist!".format(filearg))
    return filearg


def is_directory(dirarg):
    """ Type for argparse - checks that directory exists.
    """
    if not os.path.isdir(dirarg):
        raise argparse.ArgumentError(
            "The directory '{0}' does not exist!".format(dirarg))
    return dirarg


parser = argparse.ArgumentParser(description=doc)
parser.add_argument(
    "-v", "--verbose", dest="verbose", type=int, choices=[0, 1, 2], default=0,
    help="increase the verbosity level: 0 silent, [1, 2] verbose.")
parser.add_argument(
    "-c", "--subject_code", dest="subject_code", required=True,
    help="an identifier for the subject used to create the result fodler.",
    type=str)
parser.add_argument(
    "-e", "--erase", dest="erase", action="store_true",
    help="if activated, clean the result folder.")
parser.add_argument(
    "-f", "--fslconfig", dest="fslconfig", required=True, metavar="FILE",
    help="the '.sh' freesurfer configuration file.", type=is_file)
parser.add_argument(
    "-s", "--spmbin", dest="spmbin", required=True, metavar="FILE",
    help="the spm standalone file.", type=is_file)
parser.add_argument(
    "-o", "--outdir", dest="outdir", required=True, metavar="PATH",
    help="the path to output working directory.", type=is_directory)
parser.add_argument(
    "-d", "--display", dest="display", action="store_true",
    help="if activated, display the pipeline.")
parser.add_argument(
    "-n", "--normalization", dest="normalization", choices=["anat", "fmri"],
    help="slect the nomalization procedure.")
parser.add_argument(
    "-t", "--timings", dest="timings", choices=["spm", "fsl"],
    help="select the slice timing procedure.")
parser.add_argument(
    "-u", "--repetition_time", dest="repetition_time", required=True,
    help="the sequence TR.", type=float)
parser.add_argument(
    "-l", "--slice_order", dest="slice_order",  nargs='+', required=True,
    help="the sequence slice order.", type=int)
parser.add_argument(
    "-i", "--fmri", dest="fmri", required=True, metavar="FILE",
    help="the path to functional file.", type=is_file)
parser.add_argument(
    "-a", "--struct", dest="struct", metavar="FILE",
    help="the path to the structural file.", type=is_file)
parser.add_argument(
    "-r", "--template", dest="template", metavar="FILE",
    help="the path to the reference template file.", type=is_file)
parser.add_argument(
    "-k", "--template_mask", dest="template_mask", metavar="FILE",
    help="the path to the reference template mask file.", type=is_file)
parser.add_argument(
    "-w", "--ref_slice", dest="ref_slice", required=True,
    help=("the slices will be corrected to what they would have been if they "
          "were acquired when the reference slice was acquired."), type=int)
args = parser.parse_args()


"""
Check input configuration based on selected processing options.
"""
if args.normalization == "anat":
    if args.struct is None:
        raise ValueError(
            "Expect a structural file for the 'anat' normalization.")
elif args.normalization == "fmri":
    if args.template is None:
        raise ValueError(
            "Expect a template file for the 'fmri' normalization.")

"""
Create a temporary folder for capsul and a subject specific output folder.
"""
if args.verbose > 0:
    print("[info] Start SPM fMRI preprocessings...")
    print("[info] Directory: {0}.".format(args.outdir))
    print("[info] Subject: {0}.".format(args.subject_code))
    print("[info] fMRI: {0}.".format(args.fmri))
funcfile = args.fmri
sid = args.subject_code
tmp_capsul = os.path.join(args.outdir, "tmp_{0}".format(sid))
if os.path.exists(tmp_capsul):
    shutil.rmtree(tmp_capsul)
outdir = os.path.join(args.outdir, sid)

if not os.path.isdir(outdir):
    os.mkdir(outdir)
elif args.erase:
    shutil.rmtree(outdir)
    os.mkdir(outdir)
logdir = os.path.join(args.outdir, "logs")
if not os.path.isdir(logdir):
    os.mkdir(logdir)
elif args.erase:
    shutil.rmtree(logdir)
    os.mkdir(logdir)

"""
First create a study configuration.
"""
study_config = StudyConfig(
    modules=["MatlabConfig", "SPMConfig", "FSLConfig", "NipypeConfig"],
    use_smart_caching=False,
    fsl_config=args.fslconfig,
    use_fsl=True,
    use_matlab=False,
    use_spm=False,
    spm_exec=args.spmbin,
    spm_standalone=True,
    use_nipype=True,
    output_directory=tmp_capsul,
    number_of_cpus=1,
    generate_logging=True,
    use_scheduler=True)

"""
Processing definition: create the <clinfmri.preproc.FmriPreproc> that
define the different step of the processings.
"""
pipeline = get_process_instance("clinfmri.preproc.fmri_preproc.xml")

"""
It is possible to display the pipeline.
"""
if args.display:
    import sys
    from PySide import QtGui
    from capsul.qt_gui.widgets import PipelineDevelopperView

    app = QtGui.QApplication(sys.argv)
    view = PipelineDevelopperView(pipeline)
    view.show()
    app.exec_()

"""
Now to parametrize the pipeline pipeline.
"""
pipeline.fmri_file = funcfile
pipeline.realign_register_to_mean = True
pipeline.select_slicer = args.timings
pipeline.select_normalization = args.normalization
pipeline.force_repetition_time = args.repetition_time
pipeline.force_slice_orders = args.slice_order
pipeline.realign_wrap = [0, 1, 0]
pipeline.realign_write_wrap = [0, 1, 0]
pipeline.ref_slice = args.ref_slice
if args.template is not None:
    pipeline.template_file = args.template
if args.struct is not None:
    pipeline.structural_file = args.struct

"""
The pipeline is now ready to be executed.
"""
study_config.run(pipeline, executer_qc_nodes=False, verbose=1)

""" 
Organize the data with the openfmri standard.
"""
if args.normalization == "fmri" and args.timings == "spm":
    processing_struct = {
        "": [
            "clinfmri.preproc.FmriPreproc.log",
            "exitcodes_status.json",
            "parameters_status.json"],
        "bet.fsl_bet": [
            "{}_brain_mask.nii.gz",
            "{}_brain.nii.gz"],
        "spm_slicer": [
            "au{}.nii",
           "pyscript_slicetiming.m"],
        "realign": [
            "rp_au{}.txt",
            "rau{}.nii",
            "pyscript_realign.m",
            "meanau{}.nii",
            "au{}.mat"],
        "spm_normalize_template": [
            "wmeanau{}.nii",
            "pyscript_normalize.m",
            "meanau{}_sn.mat"],
        "spm_funcnormalize_template": [
            "wau{}.nii",
            "pyscript_normalize.m"],
        "smoothing": [
            "pyscript_smooth.m",
            "swau{}.nii"]
    }
else:
    raise NotImplementedError("Can't reorganize the data.")
image_basename = os.path.basename(funcfile).split(".")[0]
for step_name, step_fnames in processing_struct.items():
    for basename in step_fnames:
        basename = basename.format(image_basename)
        srcpath = os.path.join(tmp_capsul, step_name, basename)
        if srcpath.endswith(".nii"):
            gzip_file(srcpath, prefix="", output_directory=outdir)
        else:
            if ((step_name == "spm_funcnormalize_template") and
                    (basename == "pyscript_normalize.m")):
                destpath = os.path.join(outdir, "pyscript_apply_normalize.m")
            elif step_name == "":
                destpath = os.path.join(
                    logdir, "{0}_{1}".format(sid, basename))
            else:
                destpath = os.path.join(outdir, basename)
            shutil.move(srcpath, destpath)
shutil.rmtree(tmp_capsul)

"""
Create a snap to chack the registration result
"""
if args.template_mask is not None:
    snap = os.path.join(outdir, image_basename + ".pdf")
    image = os.path.join(outdir, "wmeanau{}.nii.gz".format(image_basename))
    fig = plot_mosaic(image, overlay_mask=args.template_mask,
                      title="SPM template registration")
    fig.savefig(snap, dpi=300)
