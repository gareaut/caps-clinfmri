# -*- coding: utf-8 -*-

import numpy as np
import nibabel
import Image
import os
import matplotlib.pyplot as plt
import json


def power_scores(gm_mask, before_fmri_file, after_fmri_file, displacement_file,
                 output_directory, verbose=0):
    """
    Generate QC scores and jpeg images
    from Power 2014

    CAPSUL header
    -------------
    <unit>
        <input name="gm_mask" type="File" desc="the gray matter mask to select
            relevant signal"/>
        <input name="before_fmri_file" type="File" desc="the fmri file before
            all noise corrections"/>
        <input name="after_fmri_file" type="File" desc="the fmri file after
            all noise corrections"/>
        <input name="displacement_file" type="File" desc="the json containing
            the max displacement values"/>
        <input name="output_directory" type="Directory" desc="the directory
            that will contain the jpeg file"/>
        <input name="verbose" type="Int" desc="the verbosity level"/>
        <output name="qc_image" type="File" desc="the qc image generated
            from displacement file and gm voxel timecourses"/>
    </unit>
    """

    # plot displacement
    with open(displacement_file, "r") as _file:
        mvt = json.load(_file)

    fig = plt.figure()
    ax = fig.add_subplot(3, 1, 1)
    plt.plot(mvt)
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.set_title("Max displacement")

    # load gray matter file
    mask_file = nibabel.load(gm_mask)
    mask_data = mask_file.get_data()

    if verbose > 0:
        print "{0} pixels will be tracked".format(np.sum(mask_data))

    mask_data = mask_data.ravel()

    # load data
    before_data_array = nibabel.load(before_fmri_file).get_data()

    # get voxels timecourse
    before_img = []
    for volume_nb in range(before_data_array.shape[-1]):
        vol = before_data_array[:, :, :, volume_nb].ravel()
        vol = vol[mask_data == 1]
        before_img.append(vol)

    # transform in numpy array and transpose
    before_img_array = np.asarray(before_img)
    before_img_array = np.transpose(before_img_array)

    # normalize each timecourse individually
    before_normed_img = []
    for index in range(before_img_array.shape[0]):
        line = before_img_array[index, :]
        _min = np.min(line)
        _max = np.max(line)
        line_normed = (255.0 * (line - _min) / (_max - _min)).astype(np.uint8)
        before_normed_img.append(line_normed)

    # transform back into numpy array
    before_normed_img = np.asarray(before_normed_img)

    # save image in main QC figure
    ax = fig.add_subplot(3, 1, 2)
    plt.imshow(before_normed_img, cmap="gray", aspect="auto")
    plt.setp(ax.get_xticklabels(), visible=False)
    ax.set_title("GM voxels timecourse before any correction")

    # save voxels timecourse separately
    before_img_png = Image.fromarray(before_normed_img, 'L')
    before_img_png.save(os.path.join(output_directory, "before.png"), "PNG", )

    before_img = os.path.join(output_directory, "before.png")

    # same operations with the image AFTER noise corrections
    after_data_array = nibabel.load(after_fmri_file).get_data()

    after_img = []

    for volume_nb in range(after_data_array.shape[-1]):
        vol = after_data_array[:, :, :, volume_nb].ravel()
        vol = vol[mask_data == 1]
        after_img.append(vol)

    after_img_array = np.asarray(after_img)

    after_normed_img = []
    after_img_array = np.transpose(after_img_array)

    for index in range(after_img_array.shape[0]):
        line = after_img_array[index, :]
        _min = np.min(line)
        _max = np.max(line)
        line_normed = (255.0 * (line - _min) / (_max - _min)).astype(np.uint8)
        after_normed_img.append(line_normed)

    after_normed_img = np.asarray(after_normed_img)

    # plot image
    ax = fig.add_subplot(3, 1, 3)
    plt.imshow(after_normed_img, cmap="gray", aspect="auto")
    ax.set_title("GM voxels timecourse after noise/displacement corrections")

    # save individual image
    after_img_png = Image.fromarray(after_normed_img, 'L')
    after_img_png.save(os.path.join(output_directory, "before.png"), "PNG", )
    after_img = os.path.join(output_directory, "before.png")

    # generate QC file and return
    plt.savefig(os.path.join(output_directory, "qc.png"))
    qc_image = os.path.join(output_directory, "qc.png")

    return qc_image
