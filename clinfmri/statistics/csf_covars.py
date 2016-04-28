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
import copy

# Nipy import
import nipy
from nipy.algorithms.resample import resample
from nipy.core.reference import coordinate_map as cmap
from nipy.core.reference.coordinate_system import CoordSysMaker

# Clindmri import
from clindmri.segmentation.freesurfer import mri_convert, mri_binarize

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


def resample_image(source_file, target_file, output_directory,
                   w2wmap_file=None, erode_path_nb=0, order=3,
                   cval=0.0, verbose=0):
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
    erode_path_nb: Int (optional, default 1)
        the number of path of the erosion. Performed before resampling.
    verbose: int (optional, default 0)
        the verbosity level.

    Returns
    -------
    resampled_file: str
        the resampled image.

    CAPSUL header
    -------------
    
    <process capsul_xml="2.0">
      <input name="source_file" type="file" doc="the image to resample."/>
      <input name="target_file" type="file" doc="the reference image."/>
      <input name="output_directory" type="directory" doc="the folder where the resampled image will be saved."/>
      <input name="w2wmap_file" type="file" doc="physical to physical transformation file."/>
      <input name="erode_path_nb" type="int" doc="the number of path in erosion of the mask"/>
      <input name="order" type="int" doc="interpolation mode, 0 = nearest neighbour"/>
      <input name="cval" type="float"/>
      <input name="verbose" type="int" doc="verbosity level"/>
      <output name="resampled_file" type="file" doc="the resampled image."/>
    </process>
    

    """
    # get world to world transformation
    # SPM version
    w2wmap = scipy.io.loadmat(w2wmap_file)["Affine"]
    # fsl version
#    w2wmap = np.fromfile(w2wmap_file, sep=" ")
#    w2wmap = w2wmap.reshape(4, 4)

    w2wmap = np.linalg.inv(w2wmap)

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

    # get anatomic mask
    source_image = nibabel.load(source_file)
    source_data = source_image.get_data()

    # erode
    if erode_path_nb > 0:
        eroded_image = binary_erosion(
            source_data,
            iterations=erode_path_nb).astype(source_data.dtype)

        # save
        _temp = nibabel.Nifti1Image(eroded_image, source_image.get_affine())
        eroded_file = os.path.join(output_directory,
                                   'eroded_anat_mask.nii.gz')
        nibabel.save(_temp, eroded_file)
        # Get eroded anatomic mask
        source_image = nipy.load_image(eroded_file)

    # resample
    resampled_image = resample(
        source_image, cmap_out, w2wmap, onto_shape, order=order, cval=cval)

    # save
    resampled_file = os.path.join(
        output_directory,
        "resampled_{0}".format(os.path.basename(source_file)))
    nipy.save_image(resampled_image, resampled_file)

    return resampled_file


def get_covars(csfmask_file, func_file, min_nb_of_voxels=20, nb_covars=5,
               verbose=0, output_directory=None):
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
    output_directory: str (optional, default None)
        for debuging purpose: if the verbosity level is > 1 save the mask used
        to select the functional time series.

    Returns
    -------
    covars: array (T, nb_covars)
        the requested number of covariates that represent the CSF variability.

    CAPSUL HEADER
    -------------

    
    <process capsul_xml="2.0">
      <input name="csfmask_file" type="file" doc="a binary mask of the CSF in the functional space."/>
      <input name="func_file" type="file" doc="a functional volume of size (X, Y, Z, T)."/>
      <input name="min_nb_of_voxels" type="int" doc="the criterion used to select a CSF ROI with specific size. Optional (default=50)"/>
      <input name="nb_covars" type="int" doc="the number of covariates used to explain the CSF variability. Optional (default=5)"/>
      <input name="verbose" type="int" doc="the verbosity level. Optional (default=0)"/>
      <input name="output_directory" type="directory" doc="for debuging purpose: if the verbosity level is &gt; 1 save the mask used to select the functional time series."/>
      <output name="covars" type="file" doc="the requested number of covariates that represent the CSF variability."/>
    </process>
    
    """
    # Erode the mask until we have N < min_nb_of_voxels ones
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
        csf_mask = nibabel.Nifti1Image(csf_array.astype(int),
                                       csf_image.get_affine())
        nibabel.save(csf_mask, os.path.join(output_directory,
                                            "covars_mask.nii.gz"))

    # Compute a SVD
    func_array = nibabel.load(func_file).get_data()
    csftimeseries = func_array[np.where(csf_array == 1)].T
    csftimeseries -= csftimeseries.mean(axis=0)
    u, s, v = np.linalg.svd(csftimeseries, full_matrices=False)
    if verbose > 2:
        plt.plot(s)
        plt.show()

    # Get the covariates that represent the CSF variability
    covars = u[:, :nb_covars]

    np.savetxt(os.path.join(output_directory, "covars.txt"), covars)
    covars = os.path.join(output_directory, "covars.txt")

    return covars


def complete_regressors_file(input_file, covars_file, output_directory,
                             add_extra_mvt_reg=False):
    """
    Complete the rp files with covars from an extra file

    CAPSUL HEADER
    -------------

    
    <process capsul_xml="2.0">
      <input name="input_file" type="file" doc="the current regressors file"/>
      <input name="covars_file" type="file" doc="the regressors to add file"/>
      <input name="output_directory" type="directory" doc="the directory that will contain the output file"/>
      <output name="covars" type="file" doc="the completed covars file"/>
    </process>
    
    """
    covars = np.loadtxt(covars_file)

    rp = np.load(input_file)

    if add_extra_mvt_reg:
        # Add translation square and cube, plus shift
        t, r = rp[:, :3], rp[:, 3:]
        kl = np.vstack((t[:1, :], t[:-1, :]))
        ke = np.vstack((t[1:, :], t[-1:, :]))
        out = np.column_stack((t, r, t**2, t**3, ke, kl, covars))
    else:
        out = np.column_stack((rp, covars))

    out = (out - out.mean(axis=0)) / out.std(axis=0)
    covars = os.path.join(output_directory, "complete_reg_file.txt")
    np.savetxt(covars, out, fmt="%5.8f")
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


def select_white_matter(fsdir, subjectid, mat_file, outdir, verbose=0):
    """ Select only the white matter using FreeSurfer segmentation.

    Parameters
    ----------
    fsdir: str (mandatory)
        the freesurfer working directory with all the subjects.
    subjectid: str (mandatory)
        the subject identifier.
    outdir: str (mandatory)
        the folder where the computed white matter mask image will be saved.
    verbose: int (optional, default 0)
        the verbosity level.

    Returns
    -------
    wm_file: str
        a binary segmantation of the white matter.
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

    # Get the white matter label
    wm_image = nibabel.load(wm_file)
    wm_array = wm_image.get_data()
    wm_array[np.where(wm_array != 110)] = 0
    wm_array[np.where(wm_array == 110)] = 1

    # Save the ventricules images
    wm_image = nibabel.Nifti1Image(wm_array.astype(int),
                                   wm_image.get_affine())
    wm_file = os.path.join(outdir, "white_matter.nii.gz")
    nibabel.save(wm_image, wm_file)

    # Load the affine deformation
    affine = scipy.io.loadmat(mat_file)["Affine"]
    affine = np.linalg.inv(affine)
    if verbose > 0:
        print("Affine :", affine)

    return wm_file, affine


def select_grey_matter(fsdir, subjectid, mat_file, outdir, verbose=0):
    """ Select only the grey matter using FreeSurfer segmentation.

    Parameters
    ----------
    fsdir: str (mandatory)
        the freesurfer working directory with all the subjects.
    subjectid: str (mandatory)
        the subject identifier.
    outdir: str (mandatory)
        the folder where the computed grey matter mask image will be saved.
    verbose: int (optional, default 0)
        the verbosity level.

    Returns
    -------
    wm_file: str
        a binary segmantation of the white matter.
    affine: array (4, 4)
        the physical to physical t1 to functional transformation.
    """
    # Checks
    if not os.path.isdir(outdir):
        raise ValueError("'{0}' is not a valid folder.".format(outdir))
    # Convert the segmention to the native space
    regex = os.path.join(subjectid, "mri", "aseg.mgz")
    gm_files = mri_convert(
        fsdir, regex, reslice=True, interpolation="nearest",
        fsconfig="/i2bm/local/freesurfer/SetUpFreeSurfer.sh")
    if len(gm_files) != 1:
        raise ValueError("Expect one file '{0}'.".format(gm_files))
    gm_file = gm_files[0]

    # Get the grey matter label
    gm_image = nibabel.load(gm_file)

    gm_array = gm_image.get_data()
    gm_array_tmp = copy.deepcopy(gm_array)
    gm_array[np.where(gm_array != 3)] = 0
    gm_array[np.where(gm_array == 3)] = 1
    gm_array_tmp[np.where(gm_array_tmp != 42)] = 0
    gm_array_tmp[np.where(gm_array_tmp == 42)] = 1

    gm_array += gm_array_tmp

    # Save the ventricules images
    gm_image = nibabel.Nifti1Image(gm_array.astype(int),
                                   gm_image.get_affine())
    gm_file = os.path.join(outdir, "grey_matter.nii.gz")
    nibabel.save(gm_image, gm_file)

    # Load the affine deformation
    affine = scipy.io.loadmat(mat_file)["Affine"]
    affine = np.linalg.inv(affine)
    if verbose > 0:
        print("Affine :", affine)

    return gm_file, affine


def get_freesurfer_binary_mask(fsdir, sid, outdir, outfile_name,
                               region_ids=[], regions_label=""):
    """
    Generate binary mask from freesurfer segmentation calling freesurfer
    mri_binarize function
    output a nifti file
    """

    mask = mri_binarize(fsdir, sid, outfile_name, output_directory=outdir,
                        region_ids=region_ids, regions_label=regions_label)

    mask_file = nibabel.load(mask)
    mask_file_path = os.path.join(outdir, "{0}.nii.gz".format(outfile_name))
    nibabel.save(mask_file, mask_file_path)

    return mask_file_path


def get_spm_affine_param(spm_mat_file, verbose=0):
    """
    return the spm affine parameter
    """
    affine = scipy.io.loadmat(spm_mat_file)["Affine"]
    affine = np.linalg.inv(affine)

    if verbose > 0:
        print("Affine :", affine)

    return affine


def get_spm_binary_mask(probabilistic_mask, outdir, out_mask_name,
                        threshold=0.9, erode_path=3):
    """
    Generate binary mask from spm segmentation
    """

    mask = nibabel.load(probabilistic_mask)
    mask_data = mask.get_data()
    mask_data_index = mask_data < threshold
    mask_data[mask_data_index] = 0
    mask_data_index = mask_data > 0
    mask_data[mask_data_index] = 1
    mask_data = scipy.ndimage.binary_erosion(mask_data,
                                             iterations=erode_path).astype(
        mask_data.dtype)

    mask_out = nibabel.Nifti1Image(mask_data, mask.get_affine())
    out_file = os.path.join(outdir, "{0}.nii.gz".format(out_mask_name))
    nibabel.save(mask_out, out_file)

    return out_file
