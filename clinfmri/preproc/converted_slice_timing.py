#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import os
import numpy
import nibabel


def time_serie_metadata(fmri_file, force_repetition_time=0,
                        force_slice_orders=[], slice_dim=2):
    """ Information of the time serie formatted accordingly to the selected
    software.

    <process capsul_xml="2.0">
      <input name="fmri_file" type="file" doc="the input FMRI image."/>
      <input name="force_repetition_time" type="float" doc="the sequence repetition time."/>
      <input name="force_slice_orders" type="list_int" doc="the sequence slice orders."/>
      <input name="slice_dim" type="int" doc="dimensions of a slice"/>
      <output name="repetition_time" type="float" doc="the sequence repetition time."/>
      <output name="acquisition_time" type="float" doc="the sequence acquisition time."/>
      <output name="slice_orders" type="list_int" doc="the sequence slice orders."/>
      <output name="number_of_slices" type="int" doc="the number of slices in the sequence."/>
    </process>
    """
    # Load the image and get the corresponding header
    nii = nibabel.load(fmri_file)
    header = nii.get_header()
    
    print force_slice_orders
    print force_repetition_time

    # Get image information from header if necessary
    if force_repetition_time == 0 or len(force_slice_orders) == 0:
        try:
            number_of_slices = header.get_n_slices()
        except:
            number_of_slices = nii.shape[-2]
        header.set_dim_info(slice=slice_dim)
        slice_duration = header.get_slice_duration()
        slice_times = numpy.asarray(header.get_slice_times())
    else:
        number_of_slices = nii.shape[-2]

    # Get the repetition time
    if force_repetition_time == 0:
        repetition_time = (slice_duration * number_of_slices / 1000.)  # sec
    else:
        repetition_time = force_repetition_time

    # Get the slice acquisition times
    if len(force_slice_orders) == 0:
        slice_orders = numpy.round(slice_times / slice_duration).astype(int)
        slice_orders = [index + 1 for index in slice_orders]
    else:
        slice_orders = force_slice_orders

    # Compute the acquisition time
    acquisition_time = repetition_time * (1. - 1. / float(number_of_slices))

    return repetition_time, acquisition_time, slice_orders, number_of_slices


def fsl_save_custom_timings(slice_orders, output_directory):
    """ Save acquisition slice timings in order to perform the slice timing
    with FSL.

    <process capsul_xml="2.0">
      <input name="slice_orders" type="list_int" doc="the sequence slice orders."/>
      <input name="output_directory" type="directory" doc="the output directory where fsl slice times are saved."/>
      <return name="timings_file" type="file" doc="the acquisition slice timings."/>
    </process>
    """
    # Check that the outdir is valid
    if not os.path.isdir(output_directory):
        raise ValueError("'{0}' is not a valid directory.".format(
            output_directory))

    # Type conversion
    custom_timings = (numpy.asarray(slice_orders).astype(numpy.single) - 1)

    # Express slice_times as fractions of TR
    custom_timings = custom_timings / numpy.max(custom_timings)

    # In FSL slice timing an input file is expected
    timings_file = os.path.join(output_directory, "custom_timings.txt")
    numpy.savetxt(timings_file, custom_timings)

    return timings_file
