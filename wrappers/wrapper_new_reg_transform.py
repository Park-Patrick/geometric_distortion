from nipype.interfaces.base import (
    TraitedSpec,
    CommandLine,
    CommandLineInputSpec,
    File,
    traits
)
import os

root_dir = '/data/geo_dbm/data2/'

nipype_dir = 'nipype_data/'

resampled_despot1 = 'despot1_05mm.nii.gz'
affine_xfm = '7T_to_3T_affine_xfm.txt'
affine_displacement_field_gz = 'affine_displacement_field.nii.gz'

class new_reg_transformInputSpec(CommandLineInputSpec):
    reference_img = File(argstr='-ref %s', mandatory=True,
                         position=0, desc='reference image')

    disp = File(argstr='-disp %s',
                              position = 1, desc = 'disp')
    
    disp2 = File(argstr='%s',
                            position = 2, desc = 'disp2')
                         

class new_reg_transformOutputSpec(TraitedSpec):
    output_disp2 = File(desc='output disp2')

class new_reg_transform(CommandLine):
    _cmd = 'new_reg_transform'
    input_spec = new_reg_transformInputSpec
    output_spec = new_reg_transformOutputSpec
    
    def _list_outputs(self):
            outputs = self.output_spec().get()
            
            output_disp2 = self.inputs.disp2
            outputs['output_disp2'] = output_disp2
            return outputs

if __name__ == '__main__':

    test_newtransform = new_reg_transform()
    test_newtransform.inputs.reference_img = root_dir + nipype_dir + resampled_despot1
    test_newtransform.inputs.disp = root_dir + nipype_dir + affine_xfm
    test_newtransform.inputs.disp2 = root_dir + nipype_dir + affine_displacement_field_gz
    print test_newtransform.cmdline
    #test_newtransform.run()
