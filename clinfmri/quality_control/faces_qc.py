
import json
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import nibabel as nib
import os
import numpy


def compute_scores(input_image, atlas, atlas_names, output_directory):
    """
    QC of a faces task nifti, following the atlas of face region (regions in
    wich we are supposed to find activation events)
    See Zhen 2015 NeuroImage doi:10.1016/j.neuroimage.2015.03.010
    See also the resource :
        http://www.brainactivityatlas.org/atlas/atlas-explorer/

    <unit>
        <input name="input_image" type="File" desc="the image to QC"/>
        <input name="atlas" type="List" content="File" desc="The atlas files"/>
        <input name="atlas_names" type="List" content="String"
            desc="The atlas names"/>
        <input name="output_directory" type="Directory"
            desc="The output dir"/>
        <output name="histogram" type="File"
            desc="histograms of the overlapping regions"/>
        <output name="score_file" type="File"
            desc="the json files containing all scores"/>
    </unit>

    NOTE:

    """
    # get number of iteration
    print input_image
    data_array = nib.load(input_image).get_data()
    data_array[~numpy.isnan(data_array)] = 1
    data_array[numpy.isnan(data_array)] = 0

    img = nib.Nifti1Image(data_array, nib.load(input_image).get_affine())
    nib.save(img, os.path.join(output_directory, 'debug.nii.gz'))

    scores = {}
    pp = PdfPages(os.path.join(output_directory, 'histograms.pdf'))
    for region_name, region in zip(atlas_names, atlas):

        fig = plt.figure()

        # get the atlas data
        region = nib.load(region).get_data()

        # mask
        submask_region = region > 0.1
        submask_data = data_array

        mask = numpy.multiply(submask_region.astype(int),
                              submask_data)

        # correlation
        corr = numpy.ma.masked_array(region, mask=1 - mask)
        masked_region = numpy.ma.masked_array(region, mask=1 - submask_region)

#        masked_region = numpy.ma.masked_array(region, region > 0.01)

#        to_plot = [masked_region.compressed().ravel(),
#                                               masked_corr.compressed().ravel()]))

        if numpy.sum(mask) > 0:
            plt.subplot(111)
            plt.hist([corr.compressed().ravel(),
                      masked_region.compressed().ravel()],
                     20,
                     histtype='bar',
                     color=['blue', 'red'],
                     label=["Correlation", "Atlas"]
                     )
            plt.legend()
            plt.title(region_name)

            scores["mean_{0}".format(region_name)] = corr.mean()
            scores["std_{0}".format(region_name)] = corr.std()
            scores["max_{0}".format(region_name)] = corr.max()
            scores["min_{0}".format(region_name)] = corr.min()

            pp.savefig(fig)

    pp.close()

    # write scores
    with open(os.path.join(output_directory, "scores.json"), "w") as _file:
        json.dumps(scores, _file, indent=4)

    score_file = os.path.join(output_directory, "scores.json")
    histogram = os.path.join(output_directory, 'histograms.pdf')

    return histogram, score_file
