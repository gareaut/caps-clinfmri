#! /usr/bin/env python
##########################################################################
# NSAP - Copyright (C) CEA, 2015
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
from __future__ import print_function
import os
import numpy as np
import nibabel
from scipy.ndimage.morphology import binary_erosion
import scipy.io

# Nipy import
import nipy
from nipy.algorithms.resample import resample
from nipy.core.reference import coordinate_map as cmap
from nipy.core.reference.coordinate_system import CoordSysMaker

# Clindmri import
from clindmri.segmentation.freesurfer import mri_convert

# Matplotlib import
import matplotlib.pyplot as plt


def xyz_affine(big_aff, xyz=[0, 1, 2], verbose=0):
    """ Select the xyz part of an affine transform.

    Parameters
    ----------
    big_affine: array (N, N)
        a N-1d affine matrix, N >= 4.
    xyz: list of int (optional, default [0, 1, 2])
        the x, y, z indices used to construct the desired affine matrix.
    verbose: int (optional, default 0)
        the verbosity level.

    Returns
    -------
    affine: array (4, 4)
        the desired 3d affine matrix.
    """
    # Select the matrix components of insterest
    affine = big_aff[xyz, :][:, xyz]

    # Get the associated translation
    trans = big_aff[xyz, -1].reshape((affine.shape[0], 1))

    # Contruct the final affine transformation for homogeneous coordianates
    last_row = np.zeros((1, affine.shape[1] + 1))
    last_row[0, -1] = 1.
    affine = np.hstack((affine, trans))
    affine = np.vstack((affine, last_row))
    if verbose:
        print("\nold affine:\n", big_aff, "\nnew affine:\n", affine)

    return affine


def resample_image(source_file, target_file, outdir, w2wmap=None, order=3,
                   cval=0, verbose=0):
    """ Resample the source image to match the target image using Nipy.

    Parameters
    ----------
    source_file: str (mandatory)
        the image to resample.
    target_file: str (mandatory)
        the reference image.
    outdir: str (mandatory)
        the folder where the resampled image will be saved.
    w2wmap: array (4, 4) or callable
        physical to physical transformation.
    verbose: int (optional, default 0)
        the verbosity level.

    Returns
    -------
    resampled_file: str
        the resampled image.
    """
    # Get target image information
    target_image = nipy.load_image(target_file)
    onto_shape = target_image.shape[:3]
    onto_aff = xyz_affine(target_image.affine, xyz=[0, 1, 2], verbose=verbose)

    # Define index and physical coordinate systems
    arraycoo = "ijklmnopq"[:len(onto_shape)]
    spacecoo = "xyztrsuvw"[:len(onto_shape)]
    if verbose > 0:
        print("\narraycoo: ", arraycoo, "\nspacecoo: ", spacecoo,
              "\nonto_aff\n", onto_aff)
    dmaker = CoordSysMaker(arraycoo, 'generic-array')
    rmaker = CoordSysMaker(spacecoo, 'generic-scanner')
    cm_maker = cmap.CoordMapMaker(dmaker, rmaker)
    cmap_out = cm_maker.make_affine(onto_aff)
    if verbose > 0:
        print("cmap_out:\n", cmap_out)

    # Define the default physical to physical transformation
    if w2wmap is None:
        w2wmap = np.eye(onto_aff.shape[0])
    if verbose > 0:
        print("w2wmap:\n", w2wmap)

    # Resample
    source_image = nipy.load_image(source_file)
    resampled_image = resample(
        source_image, cmap_out, w2wmap, onto_shape, order=order, cval=cval)

    # Save the resampled image
    resampled_file = os.path.join(
        outdir, "resampled_{0}".format(os.path.basename(source_file)))
    nipy.save_image(resampled_image, resampled_file)

    return resampled_file


def csf_covars(csfmask_file, func_file, min_nb_of_voxels=20, nb_covars=5,
               verbose=0, outdir=None):
    """ Compute covariates that represent the CSF variability in a functional
    timeserie.

    Parameters
    ----------
    csfmask_file: str (mandatory)
        a binary mask of the CSF in the functional space.
    func_file: str (mandatory)
        a functional volume of size (X, Y, Z, T).
    min_nb_of_voxels: int (optional, default 50)
        the criterion used to select a CSF ROI with specific size.
    nb_covars: int (optional, default 5)
        the number of covariates used to explain the CSF variability.
    verbose: int (optional, default 0)
        the verbosity level.
    outdir: str (optional, default None)
        for debuging purpose: if the verbosity level is > 1 save the mask used
        to select the functional time series.

    Returns
    -------
    covars: array (T, nb_covars)
        the requested number of covariates that represent the CSF variability.
    """
    # Erode the mask in order to have > min_nb_of_voxels ones
    csf_image = nibabel.load(csfmask_file)
    csf_array = csf_image.get_data()
    csf_array[np.where(csf_array != 0)] = 1
    if len(np.where(csf_array == 1)[0]) > min_nb_of_voxels:
        while True:
            csf_tmp_array = binary_erosion(csf_array, iterations=1)
            nb_of_ones = len(np.where(csf_tmp_array == 1)[0])
            if nb_of_ones < min_nb_of_voxels:
                break
            csf_array = csf_tmp_array
    else:
        raise ValueError(
            "Not enough CSF voxels in mask '{0}'.".format(csfmask_file))
    if verbose > 1:
        if outdir is not None:
            csf_mask = nibabel.Nifti1Image(csf_array.astype(int),
                                           csf_image.get_affine())
            nibabel.save(csf_mask, os.path.join(outdir, "covars_mask.nii.gz"))

    # Compute a SVD
    func_array = nibabel.load(func_file).get_data()
    csftimeseries = func_array[np.where(csf_array == 1)].T
    csftimeseries -= csftimeseries.mean(axis=0)
    u, s, v = np.linalg.svd(csftimeseries, full_matrices=False)
    if verbose > 1:
        plt.plot(s)
        plt.show()

    # Get the covariates that represent the CSF variability
    covars = u[:, :nb_covars]

    return covars


def select_ventricules(fsdir, subjectid, mat_file, outdir, verbose=0):
    """ Select only the ventricules using FreeSurfer segmentation.

    Parameters
    ----------
    fsdir: str (mandatory)
        the freesurfer working directory with all the subjects.
    subjectid: str (mandatory)
        the subject identifier.
    outdir: str (mandatory)
        the folder where the computed ventricules mask image will be saved.
    verbose: int (optional, default 0)
        the verbosity level.

    Returns
    -------
    ventricules_file: str
        a binary segmantation of the ventricules.
    affine: array (4, 4)
        the physical to physical t1 to functional transformation.
    """
    # Checks
    if not os.path.isdir(outdir):
        raise ValueError("'{0}' is not a valid folder.".format(outdir))

    # Convert the segmention to the native space
    regex = os.path.join(subjectid, "mri", "wm.mgz")
    wm_files = mri_convert(
        fsdir, regex, reslice=True, interpolation="nearest",
        fsconfig="/i2bm/local/freesurfer/SetUpFreeSurfer.sh")
    if len(wm_files) != 1:
        raise ValueError("Expect one file '{0}'.".format(wm_files))
    wm_file = wm_files[0]

    # Get the ventricule label
    wm_image = nibabel.load(wm_file)
    wm_array = wm_image.get_data()
    wm_array[np.where(wm_array != 250)] = 0
    wm_array[np.where(wm_array == 250)] = 1

    # Save the ventricules images
    ventricules_image = nibabel.Nifti1Image(wm_array.astype(int),
                                            wm_image.get_affine())
    ventricules_file = os.path.join(outdir, "ventricules.nii.gz")
    nibabel.save(ventricules_image, ventricules_file)

    # Load the affine deformation
    affine = scipy.io.loadmat(mat_file)["Affine"]
    affine = np.linalg.inv(affine)
    if verbose > 0:
        print("Affine :", affine)

    return ventricules_file, affine
