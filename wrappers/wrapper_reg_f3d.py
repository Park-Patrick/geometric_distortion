#nifti f3d wrapper
#reg_f3d -ref $resampled_despot1 -flo $resampled_highres_bet -rigOnly -aff $rigid_xfm -res $rigid_nii
from nipype.interfaces.base import(
    TraitedSpec,
    CommandLineInputSpec,
    CommandLine,
    File,
    traits
)
import os

#TODO: Create text file to import file names from
root_dir = '/home/ROBARTS/ppark/Documents/park_test_geo_dbm/nipype/EPI_P021/'

resampled_despot1 = 'despot1_05mm.nii.gz'
resampled_highres_bet = 'MPRAGE_std_N4_bet_05mm.nii.gz'

rigid_xfm = '7T_to_3T_rigid_xfm.txt'
rigid_nii = '7T_to_3T_rigid.nii'

affine_xfm = '7T_to_3T_affine_xfm.txt'
affine_nii = '7T_to_3T_affine.nii'

nlin_cpp = '7T_to_3T_nlin_cpp.nii'
nlin_nii = '7T_to_3T_nlin.nii'

class reg_f3dInputSpec(CommandLineInputSpec):
    reference_img = File(argstr='-ref %s', mandatory=True,
                         position=0, desc='reference image')
    floating_img = File(argstr='-flo %s', mandatory=True,
                        position=1, desc='floating image')
    option_affine = File(argstr='-aff %s',
                         position = 2, desc = 'output affine transformation')
    control_point_grid = File(argstr='-cpp %s',
                              position = 3, desc = 'control_point_grid filename')
    resampled_img = File(argstr='-res %s',
                         position = 4, desc = 'resampled image')
    
                         

class reg_f3dOutputSpec(TraitedSpec):
    outfile_cpp = File(exists=True, desc='outfile_cpp')
    outfile_nii = File(exists=True, desc='outfile_nii')

class reg_f3d(CommandLine):
    _cmd = 'reg_f3d'
    input_spec = reg_f3dInputSpec
    output_spec = reg_f3dOutputSpec
    
    def _list_outputs(self):
            outputs = self.output_spec().get()
            outfile_cpp = self.inputs.control_point_grid
            outputs['outfile_cpp'] = outfile_cpp
            
            outfile_nii = self.inputs.resampled_img
            outputs['outfile_nii'] = outfile_nii
#            outputs['out_file'] = os.path.abspath(self.inputs.input_file + "_registered")
            return outputs

if __name__ == '__main__':

    test_f3d = reg_f3d()
    test_f3d.inputs.reference_img = root_dir + resampled_despot1
    test_f3d.inputs.floating_img = root_dir + resampled_highres_bet
    test_f3d.inputs.option_affine = root_dir + affine_xfm
    test_f3d.inputs.control_point_grid = root_dir + nlin_cpp
    test_f3d.inputs.resampled_img = root_dir + nlin_nii
    print test_f3d.cmdline
#    test_f3d.run()
