#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
from __future__ import division
from __future__ import print_function
import os
import json
import nibabel
import numpy

# ITK import
import SimpleITK as sitk

# MONKEYFMRI import
from pclinfmri.utils.runtime import progress_bar


def command_iteration(method):
    """ ITK command used to dispalay the optimization status.
    """
    print("{0:3} = {1:7.5f} : {2}".format(method.GetOptimizerIteration(),
                                          method.GetMetricValue(),
                                          method.GetOptimizerPosition()))


def slice_registration(in_file, slice_shift=5, registration_prefix="w",
                       transformation_prefix="rp", output_directory=None,
                       verbose=0):
    """ Register slices in an EPI 4D volume.

    Compute the meadian image and regsiter the 4d volume image slices
    to the median image slices.

    Parameters
    ----------
    in_file: str (mandatory)
        the input 4D EPI volume.
    slice_shift: int (optional, default 5)
        slice to exclude (+slice_shift, dim - slice_shift).
    registration_prefix: str (optional, default 'w')
        the prefix of the output registered volume file.
    transformation_prefix: str (optional, default 'rp')
        the prefix of the output parameters file.
    output_directory: str (optional, default None)
        the output directory where the rectified image is saved.
        If None use the same directory as the input image.
    verbose: int (optional, default 0)
        parameter to control the verbosity of the code.

    Returns
    -------
    register_file: str
        the output temporal slice registered volume.
    transfromation_file: str
        the computed slice to slice 2d affine transformation parameters.

    Example
    -------

    >>> import pclinfmri.utils.registration import slice_registration
    >>> reg = slice_registration('image.nii', 5, 'w', 'rp', None, 1)

    <process>
        <return name="register_file" type="File" desc="the output temporal
            slice registered volume."/>
        <return name="transformation_file" type="File" desc="the computed
            slice to slice 2d affine transformation parameters."/>
        <input name="in_file" type="File" desc="the input 4D EPI volume."/>
        <input name="slice_shift" type="Int" desc="slice to exclude
            (+slice_shift, dim - slice_shift)."/>
        <input name="registration_prefix" type="String" desc="the prefix of
            the output registered volume file."/>
        <input name="transformation_prefix" type="String" desc="the prefix of
            the output parameters file."/>
        <input name="output_directory" type="Directory" desc="the output
            directory where the rectified image is saved."/>
    </process>
    """
    # Check the input image exists on the file system
    if not os.path.isfile(in_file):
        raise ValueError("'{0}' is not a valid filename.".format(in_file))

    # Check that the outdir is valid
    if output_directory is not None:
        if not os.path.isdir(output_directory):
            raise ValueError("'{0}' is not a valid directory.".format(
                output_directory))
    else:
        output_directory = os.path.dirname(in_file)

    # Load the 4d EPI volume and get the associated data
    time_serie_image = nibabel.load(in_file)
    time_serie_data = time_serie_image.get_data()
    time_serie_volumes = time_serie_data.shape[-1]
    time_serie_slices = time_serie_data.shape[-2]

    # Build the median image from the 4d volume
    median_data = numpy.median(time_serie_data, 3).astype(numpy.single)

    # Set up the registration
    registration = sitk.ImageRegistrationMethod()
    registration.SetInterpolator(sitk.sitkLinear)
    registration.SetMetricAsMeanSquares()
    registration.SetOptimizerAsRegularStepGradientDescent(
        learningRate=1.0, minStep=.001, numberOfIterations=300,
        relaxationFactor=0.8)
    registration.SetOptimizerScalesFromIndexShift()
    if verbose > 1:
        registration.AddCommand(sitk.sitkIterationEvent,
                                lambda: command_iteration(registration))

    # Go through all volumes
    reg_parameters = {}
    ts_reg_data = numpy.zeros(time_serie_data.shape, dtype=numpy.single)
    for vol_index in range(time_serie_volumes):

        # Display a progress bar dependending on the verbisity level
        if verbose > 0:
            progress_bar((vol_index + 1) / time_serie_volumes,
                         "Volume '{0}' Registration".format(vol_index))

        # Go through all slices
        for slice_index in range(slice_shift, time_serie_slices - slice_shift):

            # Get the current slice array data
            median_slice_data = median_data[..., slice_index]
            ts_slice_data = time_serie_data[..., slice_index, vol_index]

            # Build itk image from the previous slice array data
            itk_fixed_image = sitk.GetImageFromArray(median_slice_data)
            itk_moving_image = sitk.GetImageFromArray(ts_slice_data)
            itk_fixed_image = sitk.Normalize(itk_fixed_image)
            itk_moving_image = sitk.Normalize(itk_moving_image)

            # Set initial transform
            tx = sitk.Transform(2, sitk.sitkAffine)
            tx.SetIdentity()
            registration.SetInitialTransform(tx)

            # Execute the registration
            final_transform = registration.Execute(itk_fixed_image,
                                                   itk_moving_image)
            if verbose > 1:
                print("-------")
                print(final_transform)
                print("Optimizer stop condition: {0}".format(
                    registration.GetOptimizerStopConditionDescription()))
                print(" Iteration: {0}".format(
                    registration.GetOptimizerIteration()))
                print(" Metric value: {0}".format(
                    registration.GetMetricValue()))

            # Resample the time serie image
            resampler = sitk.ResampleImageFilter()
            resampler.SetReferenceImage(itk_fixed_image)
            resampler.SetInterpolator(sitk.sitkBSpline)
            resampler.SetDefaultPixelValue(0)
            resampler.SetTransform(final_transform)
            itk_resample_image = resampler.Execute(itk_moving_image)
            itk_resample_image = sitk.Cast(
                sitk.RescaleIntensity(itk_resample_image), sitk.sitkUInt16)

            # Get the data array of the resmpled image
            ts_reg_data[..., slice_index, vol_index] = (
                sitk.GetArrayFromImage(itk_resample_image))

            # Save the transform parameters
            reg_parameters["{0}:{1}".format(vol_index, slice_index)] = (
                final_transform.GetParameters())

    # Save the registered image and associated transformations
    fsplit = os.path.split(in_file)
    register_file = os.path.join(output_directory, registration_prefix + 
                                 fsplit[1])
    transfromation_file = os.path.join(output_directory, transformation_prefix +
                                       fsplit[1].split(".")[0] + ".json")
    time_serie_image = nibabel.Nifti1Image(
        ts_reg_data, time_serie_image.get_affine())
    nibabel.save(time_serie_image, register_file)
    with open(transfromation_file, "wb") as openfile:
        json.dump(reg_parameters, openfile, sort_keys=True, indent=4)

    return register_file, transfromation_file
