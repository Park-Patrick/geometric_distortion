#nifti aladin wrapper
#reg_aladin -ref $resampled_despot1 -flo $resampled_highres_bet -rigOnly -aff $rigid_xfm -res $rigid_nii
from nipype.interfaces.base import(
    TraitedSpec,
    CommandLineInputSpec,
    CommandLine,
    File,
    traits
)
import os

root_dir = '/home/ROBARTS/ppark/Documents/park_test_geo_dbm/nipype/EPI_P021/'

resampled_despot1 = 'despot1_05mm.nii.gz'
resampled_highres_bet = 'MPRAGE_std_N4_bet_05mm.nii.gz'

rigid_xfm = '7T_to_3T_rigid_xfm.txt'
rigid_nii = '7T_to_3T_rigid.nii'

class reg_aladinInputSpec(CommandLineInputSpec):
    reference_img = File(argstr='-ref %s', mandatory=True,
                         position=0, desc='reference image')
    floating_img = File(argstr='-flo %s', mandatory=True,
                        position=1, desc='floating image')
    
    #choose one
    option_rigOnly = traits.Bool(argstr = '-rigOnly', 
                               position = 2, desc = 'rigid only flag')
    input_affine_trans = File(argstr='-inaff %s',
                         position = 2, desc = 'input affine transformation')

    option_affine = File(argstr='-aff %s',
                         position = 3, desc = 'output affine transformation')
    resampled_image = File(argstr='-res %s',
                         position = 4, desc = 'resampled image')
    

class reg_aladinOutputSpec(TraitedSpec):
    output_affine_transformation = File(exists=True, desc='output affine transformation')
    output_resampled_image = File(exists=True, desc='output resampled image')
    

class reg_aladin(CommandLine):
    _cmd = 'reg_aladin'
    input_spec = reg_aladinInputSpec
    output_spec = reg_aladinOutputSpec
    
    def _list_outputs(self):
            outputs = self.output_spec().get()
            
            #doesn't have to be two lines in this case
            #but having temp variable is useful in cases where outfile name needs to be automatically generated
            output_affine_transformation = self.inputs.option_affine
            outputs['output_affine_transformation'] = output_affine_transformation
            
            output_resampled_image = self.inputs.resampled_image
            outputs['output_resampled_image'] = output_resampled_image
#            outputs['output_affine_transform'] = os.path.abspath(self.inputs.input_file + "_registered")
#            outputs['output_resampled_image'] = os.path.join()            
            return outputs

if __name__ == '__main__':

    test_aladin = reg_aladin()
    test_aladin.inputs.reference_img = root_dir + resampled_despot1
    test_aladin.inputs.floating_img = root_dir + resampled_highres_bet
    test_aladin.inputs.option_rigOnly = True
    test_aladin.inputs.option_affine = root_dir + rigid_xfm
    test_aladin.inputs.resampled_image = root_dir + rigid_nii
    print test_aladin.cmdline
#    test_aladin.run()
    


#class FlirtInputSpec(FSLCommandInputSpec):
#    in_file = File(exists=True, argstr='-in %s', mandatory=True,
#                   position=0, desc='input file')
#    reference = File(exists=True, argstr='-ref %s', mandatory=True,
#                     position=1, desc='reference file')
#    out_file = File(argstr='-out %s', desc='registered output file',
#                    name_source=['in_file'], name_template='%s_flirt',
#                    position=2, hash_files=False)
#    out_matrix_file = File(argstr='-omat %s',
#                           name_source=['in_file'], keep_extension=True,
#                           name_template='%s_flirt.mat',
#                           desc='output affine matrix in 4x4 asciii format',
#                           position=3, hash_files=False)
#    out_log = File(name_source=['in_file'], keep_extension=True,
#                   requires=['save_log'],
#                   name_template='%s_flirt.log', desc='output log')
#
#class FlirtOutputSpec(TraitedSpec):
#    out_file = File(exists=True,
#                    desc='path/name of registered file (if generated)')
#    out_matrix_file = File(exists=True,
#                           desc='path/name of calculated affine transform '
#                           '(if generated)')
#    out_log = File(desc='path/name of output log (if generated)')
#
#class Flirt(FSLCommand):
#    _cmd = 'flirt'
#    input_spec = FlirtInputSpec
#    output_spec = FlirtOutputSpec
                    
