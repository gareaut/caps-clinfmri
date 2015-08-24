#! /usr/bin/env python
##########################################################################
# NSAp - Copyright (C) CEA, 2013
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################


"""SPM wrappers for preprocessing data
"""

__docformat__ = 'restructuredtext'

# Standard library imports
from copy import deepcopy
import os

# Third-party imports
import numpy as np

# Local imports
from nipype.interfaces.base import (OutputMultiPath, TraitedSpec, isdefined,
                                    traits, InputMultiPath, File)
from nipype.interfaces.spm.base import (SPMCommand, scans_for_fname,
                                        func_is_3d,
                                        scans_for_fnames, SPMCommandInputSpec)
from nipype.utils.filemanip import (fname_presuffix, filename_to_list,
                                    list_to_filename, split_filename)


class ApplyDeformationFieldInputSpec(SPMCommandInputSpec):
    """
    Parameters
    ----------
    in_files : list of str (mandatory)
        Files on which the deformation is applied.

    deformation_field : str (mandatory)
        SN SPM deformation file.

    bounding_box : list of list of float
        3x2-element list of lists (opt).

    voxel_sizes : list of float
        3-element list (opt).

    interpolation : int
        Degree of b-spline used for interpolation (from 0 to 7).
    """
    in_files = InputMultiPath(
        File(exists=True),
        mandatory=True,
        field='fnames',
        desc='Files on which deformation is applied')
    deformation_field = File(
        exists=True,
        mandatory=True,
        field='comp{1}.def',
        desc='SN SPM deformation file')
    bounding_box = traits.List(
        traits.List(traits.Float(),
        minlen=3, maxlen=3),
        field='comp{2}.idbbvox.bb',
        minlen=2, maxlen=2,
        desc='3x2-element list of lists (opt)')
    voxel_sizes = traits.List(
        traits.Float(),
        [1., 1., 1.],
        field='comp{2}.idbbvox.vox',
        minlen=3, maxlen=3,
        desc='3-element list (opt)')
    interpolation = traits.Range(
        low=0,
        high=7,
        field='interp',
        desc='degree of b-spline used for interpolation')
    out_prefix = traits.String(
        'w', field='prefix',
        usedefault=True,
        desc='aplly deformation field output prefix')


class ApplyDeformationFieldOutputSpec(TraitedSpec):
    """
    Returns
    -------
    normalized_files : list of str
        Converted files.
    """
    normalized_files = OutputMultiPath(
        File(exists=True),
        desc='converted files')


class ApplyDeformationField(SPMCommand):
    """ Uses SPM to apply inverse deformation field to given files.

    Examples
    --------
    >>> import nsap.nipype.spm_interfaces as spm
    >>> f = spm.ApplyformationField()
    >>> f.inputs.in_files = 'functional.nii'
    >>> f.inputs.deformation_field = 'y_t1_localizer.nii'
    >>> f.run()
    """

    input_spec = ApplyDeformationFieldInputSpec
    output_spec = ApplyDeformationFieldOutputSpec

    _jobtype = 'util'
    _jobname = 'defs'

    def _format_arg(self, opt, spec, val):
        """Convert input to appropriate format for SPM
        """
        if opt in ['in_files']:
            return scans_for_fnames(filename_to_list(val), keep4d=False)
        if opt == 'deformation_field':
            return np.array([list_to_filename(val)], dtype=object)
        return super(ApplyDeformationField, self)._format_arg(opt, spec, val)

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['normalized_files'] = []
        for filename in self.inputs.in_files:
            _, fname = os.path.split(filename)
            outputs['normalized_files'].append(
                os.path.join(self.inputs.output_directory,
                             '%s%s' % (self.inputs.out_prefix, fname)))
        return outputs



class NormalizeInputSpec(SPMCommandInputSpec):
    template = File(exists=True, field='eoptions.template',
                    desc='template file to normalize to',
                    mandatory=True, xor=['parameter_file'],
                    copyfile=False)
    source = File(exists=True, field='subj.source',
                  desc='file to normalize to template',
                  xor=['parameter_file'],
                  mandatory=True, copyfile=True)
    jobtype = traits.Enum('estwrite', 'est', 'write',
                          desc='one of: est, write, estwrite (opt, estwrite)',
                          usedefault=True)
    apply_to_files = InputMultiPath(traits.Either(File(exists=True),
                                                  traits.List(File(exists=True))),
                                    field='subj.resample',
                                    desc='files to apply transformation to (opt)',
                                    copyfile=True)
    parameter_file = File(field='subj.matname', mandatory=True,
                          xor=['source', 'template'],
                          desc='normalization parameter file*_sn.mat', copyfile=False)
    source_weight = File(field='subj.wtsrc',
                                 desc='name of weighting image for source (opt)', copyfile=False)
    template_weight = File(field='eoptions.weight',
                                   desc='name of weighting image for template (opt)', copyfile=False)
    source_image_smoothing = traits.Float(field='eoptions.smosrc',
                                          desc='source smoothing (opt)')
    template_image_smoothing = traits.Float(field='eoptions.smoref',
                                            desc='template smoothing (opt)')
    affine_regularization_type = traits.Enum('mni', 'size', 'none', field='eoptions.regype',
                                              desc='mni, size, none (opt)')
    DCT_period_cutoff = traits.Float(field='eoptions.cutoff',
                                     desc='Cutoff of for DCT bases (opt)')
    nonlinear_iterations = traits.Int(field='eoptions.nits',
                     desc='Number of iterations of nonlinear warping (opt)')
    nonlinear_regularization = traits.Float(field='eoptions.reg',
                                            desc='the amount of the regularization for the nonlinear part of the normalization (opt)')
    write_preserve = traits.Bool(field='roptions.preserve',
                     desc='True/False warped images are modulated (opt,)')
    write_bounding_box = traits.List(traits.List(traits.Float(), minlen=3,
                                                 maxlen=3),
                                     field='roptions.bb', minlen=2, maxlen=2,
                                     desc='3x2-element list of lists (opt)')
    write_voxel_sizes = traits.List(traits.Float(), field='roptions.vox',
                                    minlen=3, maxlen=3,
                                    desc='3-element list (opt)')
    write_interp = traits.Range(low=0, high=7, field='roptions.interp',
                        desc='degree of b-spline used for interpolation')
    write_wrap = traits.List(traits.Int(), field='roptions.wrap',
                        desc=('Check if interpolation should wrap in [x,y,z] '
                              '- list of bools (opt)'))
    out_prefix = traits.String('w', field='roptions.prefix', usedefault=True,
                               desc='normalized output prefix')


class NormalizeOutputSpec(TraitedSpec):
    normalization_parameters = File(exists=True, desc='MAT files containing the normalization parameters')
    normalized_source = File(exists=True, desc='Normalized source files')
    normalized_files = OutputMultiPath(File(exists=True), desc='Normalized other files')


class Normalize(SPMCommand):
    """use spm_normalise for warping an image to a template
    http://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf#page=51
    Examples
    --------
    >>> import nipype.interfaces.spm as spm
    >>> norm = spm.Normalize()
    >>> norm.inputs.source = 'functional.nii'
    >>> norm.run() # doctest: +SKIP
    """

    input_spec = NormalizeInputSpec
    output_spec = NormalizeOutputSpec
    _jobtype = 'spatial'
    _jobname = 'normalise'

    def _format_arg(self, opt, spec, val):
        """Convert input to appropriate format for spm
        """
        if opt == 'template':
            return scans_for_fname(filename_to_list(val))
        if opt == 'source':
            return scans_for_fname(filename_to_list(val))
        if opt == 'apply_to_files':
            return scans_for_fnames(filename_to_list(val))
        if opt == 'parameter_file':
            return np.array([list_to_filename(val)], dtype=object)
        if opt in ['write_wrap']:
            if len(val) != 3:
                raise ValueError('%s must have 3 elements' % opt)
        return super(Normalize, self)._format_arg(opt, spec, val)

    def _parse_inputs(self):
        """validate spm realign options if set to None ignore
        """
        einputs = super(Normalize, self)._parse_inputs(skip=('jobtype',
                                                             'apply_to_files'))
        if isdefined(self.inputs.apply_to_files):
            inputfiles = deepcopy(self.inputs.apply_to_files)
            if isdefined(self.inputs.source):
                inputfiles.extend(self.inputs.source)
            einputs[0]['subj']['resample'] = scans_for_fnames(inputfiles)
        jobtype = self.inputs.jobtype
        if jobtype in ['estwrite', 'write']:
            if not isdefined(self.inputs.apply_to_files):
                if isdefined(self.inputs.source):
                    einputs[0]['subj']['resample'] = scans_for_fname(self.inputs.source)
        return [{'%s' % (jobtype): einputs[0]}]

    def _list_outputs(self):
        outputs = self._outputs().get()

        jobtype = self.inputs.jobtype
        if jobtype.startswith('est'):
            outputs['normalization_parameters'] = []
            for imgf in filename_to_list(self.inputs.source):
                outputs['normalization_parameters'].append(fname_presuffix(imgf, suffix='_sn.mat', use_ext=False))
            outputs['normalization_parameters'] = list_to_filename(outputs['normalization_parameters'])

        if self.inputs.jobtype == "estimate":
            if isdefined(self.inputs.apply_to_files):
                outputs['normalized_files'] = self.inputs.apply_to_files
            outputs['normalized_source'] = self.inputs.source
        elif 'write' in self.inputs.jobtype:
            outputs['normalized_files'] = []
            if isdefined(self.inputs.apply_to_files):
                filelist = filename_to_list(self.inputs.apply_to_files)
                for f in filelist:
                    if isinstance(f, list):
                        run = [fname_presuffix(in_f, prefix=self.inputs.out_prefix) for in_f in f]
                    else:
                        run = [fname_presuffix(f, prefix=self.inputs.out_prefix)]
                    outputs['normalized_files'].extend(run)
            if isdefined(self.inputs.source):
                outputs['normalized_source'] = fname_presuffix(self.inputs.source, prefix=self.inputs.out_prefix)

        return outputs


class ResliceToReferenceInput(SPMCommandInputSpec):
    in_files = InputMultiPath(
        File(exists=True), mandatory=True, field='fnames',
        desc='Files on which deformation is applied')
    target = File(
        exists=True,
        field='comp{1}.id.space',
        desc='File defining target space')
    interpolation = traits.Range(
        low=0, high=7, field='interp',
        desc='degree of b-spline used for interpolation')

    bounding_box = traits.List(traits.List(traits.Float(), minlen=3, maxlen=3),
        field='comp{1}.idbbvox.bb', minlen=2, maxlen=2,
        desc='3x2-element list of lists (opt)')
    voxel_sizes = traits.List(
        traits.Float(),
        field='comp{1}.idbbvox.vox',
        minlen=3, maxlen=3,
        desc='3-element list (opt)')


class ResliceToReferenceOutput(TraitedSpec):
    out_files = OutputMultiPath(File(exists=True),
                                desc='Transformed files')


class ResliceToReference(SPMCommand):
    """ Uses spm to reslice a volume to a target image space or to a provided voxel size and bounding box
    Examples
    --------
    >>> import nipype.interfaces.spm.utils as spmu
    >>> r2ref = spmu.ResliceToReference()
    >>> r2ref.inputs.in_files = 'functional.nii'
    >>> r2ref.inputs.target = 'structural.nii'
    >>> r2ref.run() # doctest: +SKIP
    """

    input_spec = ResliceToReferenceInput
    output_spec = ResliceToReferenceOutput

    _jobtype = 'util'
    _jobname = 'defs'

    def _format_arg(self, opt, spec, val):
        """Convert input to appropriate format for spm
        """
        if opt == 'in_files':
            return scans_for_fnames(filename_to_list(val))
        if opt == 'target':
            return scans_for_fname(filename_to_list(val))
        if opt == 'deformation':
            return np.array([list_to_filename(val)], dtype=object)
        if opt == 'deformation_field':
            return np.array([list_to_filename(val)], dtype=object)
        return val

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['out_files'] = []
        for filename in self.inputs.in_files:
            _, fname = os.path.split(filename)
            outputs['out_files'].append(os.path.realpath('w%s' % fname))
        return outputs
