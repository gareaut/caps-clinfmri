#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import sys
import os
import numpy
import nibabel


def time_serie_metadata(fmri_file, force_repetition_time=None,
                        force_slice_orders=None):
    """ Information of the time serie formatted accordingly to the selected
    software.

    <process>
        <return name="repetition_time" type="Float"
            desc="the sequence repetition time."/>
        <return name="acquisition_time" type="Float"
            desc="the sequence acquisition time."/>
        <return name="slice_orders" type="List_Int"
            desc="the sequence slice orders."/>
        <return name="number_of_slices" type="Int"
            desc="the number of slices in the sequence."/>
        <input name="fmri_file" type="File" desc="the input FMRI image."/>
        <input name="force_repetition_time" type="Float" desc="the sequence
            repetition time."/>
        <input name="force_slice_orders" type="List_Int" desc="the sequence
            slice orders."/>
    </process>
    """
    # Load the image and get the corresponding header
    nii = nibabel.load(fmri_file)
    header = nii.get_header()

    # Get image information from header if necessary
    if not force_repetition_time or not force_slice_orders:
        number_of_slices = header.get_n_slices()
        slice_duration = header.get_slice_duration()
        slice_times = numpy.asarray(header.get_slice_times())
    else:
        number_of_slices = nii.shape[-2]

    # Get the repetition time
    if not force_repetition_time:
        repetition_time = (slice_duration * number_of_slices / 1000.)  # sec
    else:
        repetition_time = force_repetition_time

    # Get the slice acquisition times
    if not force_slice_orders:
        slice_orders = numpy.round(slice_times / slice_duration).astype(int)
        slice_orders = [index + 1 for index in slice_orders]
    else:
        slice_orders = force_slice_orders

    # Compute the acquisition time
    acquisition_time = repetition_time * (1. - 1. / 40.)

    return repetition_time, acquisition_time, slice_orders, number_of_slices


def fsl_save_custom_timings(slice_orders, output_directory):
    """ Save acquisition slice timings in order to perform the slice timing
    with FSL.

    <process>
        <return name="timings_file" type="File" desc="the acquisition slice
            timings."/>
        <input name="slice_orders" type="List_Int" desc="the sequence
            slice orders."/>
        <input name="output_directory" type="Directory" desc="the output
            directory where fsl slice times are saved."/>
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







