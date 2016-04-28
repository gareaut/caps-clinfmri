##########################################################################
# NSAp - Copyright (C) CEA, 2013-2015
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import os
import glob
import nibabel
import numpy
import json
import csv
import shutil
from nibabel import freesurfer

# Clindmri import
from clindmri.extensions.freesurfer.exceptions import FreeSurferRuntimeError
from clindmri.extensions.freesurfer.wrappers import FSWrapper



class MutualExclusiveArg(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class MandatoryArgMissing(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


def mri_binarize(fsdir, sid, output_directory=None,
                 interpolation="interpolate",
                 region_ids=[], regions_label='-',
                 fsconfig="/i2bm/local/freesurfer/SetUpFreeSurfer.sh"):
    """ Generate binary mask from freesurfer segmentation.
        Regions can be defined either by their name (see list below) or their
        identification number
        wm
            Set match vals to 2, 41, 77, 251-255 (aseg for cerebral WM)
        ventricles
            Set match vals those for aseg ventricles+choroid (not 4th)
        wm+vcsf
            WM and ventricular CSF, including choroid (not 4th)
        subcort-gm
            Subcortical gray matter
        then the mask is converted in nifti format and reshaped to match the
        original file's dimensions.
    
    <process capsul_xml="2.0">
      <input name="fsdir" type="directory" doc="The freesurfer working directory with all the subjects."/>
      <input name="sid" type="string" doc="The current subject identifier."/>
      <input name="output_directory" type="directory" doc="The freesurfer runtime folder."/>
      <input name="region_ids" type="list_int" doc="The identification numbers of the regions to segment"/>
      <input name="regions_label" type="string" doc="The region name, must be recognized by freesurfer"/>
      <input name="interpolation" type="string" doc="Define the interpolation method, must be from 'interpolate', 'weighted', 'nearest', 'cubic' (default='interpolate')"/>
      <input name="fsconfig" type="file" doc="The freesurfer configuration batch."/>
      <output name="mask_file" type="file" doc="The resulting mask"/>
    </process>
    
    """
    # Create the fs output directory if necessary
    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)

    # Check the interpolation method
    if interpolation not in ["interpolate", "weighted", "nearest", "cubic"]:
        raise ValueError(
            "'{0}' is not a valid interpolation method.".format(interpolation))

    # mutual exclusion check
    if len(region_ids) > 0 and regions_label != '-':
        raise MutualExclusiveArg("Both 'region_ids' and 'regions_label' are"
                                 " defined, only one of them should be")
    if not region_ids and not regions_label:
        raise MandatoryArgMissing("At least one of 'region_ids' or "
                                  "'regions_label' must be defined")

    fsfile = os.path.join(fsdir, sid, "mri", "aseg.mgz")
    reference_file = os.path.join(fsdir, sid, "mri", "rawavg.mgz")
    if not os.path.isfile(reference_file):
        raise ValueError("'{0}' does not exists, can't reslice image "
                         "'{1}'.".format(reference_file, fsfile))

    # Calls freesurfer: first generate binary mask then convert it in nifti
    # MASK generation
    # common part of the command
    if region_ids:
        mask_file = "freesurfer_{0}.mgz".format(
            "_".join([str(x) for x in region_ids]))
        out_mgz_path = os.path.join(output_directory, mask_file)
        cmd = ["mri_binarize", "--i", fsfile, "--o", out_mgz_path]
        for ids in region_ids:
            cmd.append("--match")
            cmd.append(str(ids))
    else:
        mask_file = "freesurfer_{0}.mgz".format(regions_label)
        out_mgz_path = os.path.join(output_directory, mask_file)
        regions = "--{0}".format(regions_label)
        cmd = ["mri_binarize", "--i", fsfile, "--o", out_mgz_path, regions]

    binarize = FSWrapper(cmd, shfile=fsconfig)
    binarize()
    if binarize.exitcode != 0:
        raise FreeSurferRuntimeError(
            binarize.cmd[0], " ".join(binarize.cmd[1:]),
            binarize.stderr + binarize.stdout)

    # CONVERSION
    mask_file = out_mgz_path.replace(".mgz", ".nii.gz")
    cmd = ["mri_convert", "--resample_type", interpolation,
           "--reslice_like", reference_file,
           out_mgz_path, mask_file]

    convert = FSWrapper(cmd, shfile=fsconfig)
    convert()
    if convert.exitcode != 0:
        raise FreeSurferRuntimeError(
            convert.cmd[0], " ".join(convert.cmd[1:]),
            convert.stderr + convert.stdout)

    return mask_file
