


from pclinfmri.utils.reorientation import reorient_image

in_file = "/volatile/nsap/pclinfmri/reorientation/t1.nii"
rectified_image = reorient_image(in_file, 'RIA', 's', None)
