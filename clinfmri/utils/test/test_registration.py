

from pclinfmri.utils.registration import slice_registration

in_file = "/volatile/nsap/pclinfmri/reorientation/ts.nii"
register_file, transfromation_file = slice_registration(
    in_file, slice_shift=5, registration_prefix="w",
    transformation_prefix="rp", outdir=None, verbose=1)
