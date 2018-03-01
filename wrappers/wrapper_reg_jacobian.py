from nipype.interfaces.base import (
    TraitedSpec,
    CommandLine,
    CommandLineInputSpec,
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

log_jacobian = '7T_to_3T_log_jacobian.nii'
mat_jacobian = '7T_to_3T_mat_jacobian.nii'
jacobian = '7T_to_3T_jacobian.nii'

class reg_jacobianInputSpec(CommandLineInputSpec):
    reference_img = File(argstr='-ref %s', mandatory=True,
                         position=0, desc='reference image')
    control_point_grid = File(argstr='-cpp %s',
                              position = 1, desc = 'control_point_grid filename')
    option_affine = File(argstr='-aff %s',
                         position = 2, desc = 'output affine transformation')
    jacobianLog = File(argstr='-jacL %s',
                       position = 3, desc = 'jacobianLog')
    jacobianMatrix = File(argstr='-jacM %s',
                          position = 3, desc = 'jacobianMatrix')
    jacobianDeterminant = File(argstr='-jac %s',
                               position = 3, desc = 'jacobianDeterminant')
    
                         

class reg_jacobianOutputSpec(TraitedSpec):
    outfile_jacobianLog = File(desc='outfile_jacobianLog')
    outfile_jacobianMatrix = File(desc='outfile_jacobianMatrix')
    outfile_jacobianDeterminant = File(desc='outfile_jacobianDeterminant')
    

class reg_jacobian(CommandLine):
    _cmd = 'reg_jacobian'
    input_spec = reg_jacobianInputSpec
    output_spec = reg_jacobianOutputSpec
    
    def _list_outputs(self):
            outputs = self.output_spec().get()
            outfile_jacobianLog = self.inputs.jacobianLog
            outputs['outfile_jacobianLog'] = outfile_jacobianLog
            
            outfile_jacobianMatrix = self.inputs.jacobianMatrix
            outputs['outfile_jacobianMatrix'] = outfile_jacobianMatrix
            
            outfile_jacobianDeterminant = self.inputs.jacobianDeterminant
            outputs['outfile_jacobianDeterminant'] = outfile_jacobianDeterminant
#            outputs['out_file'] = os.path.abspath(self.inputs.input_file + "_registered")
            return outputs

if __name__ == '__main__':

    test_jacobian = reg_jacobian()
    test_jacobian.inputs.reference_img = root_dir + resampled_despot1
    test_jacobian.inputs.control_point_grid = root_dir + nlin_cpp
    test_jacobian.inputs.option_affine = root_dir + affine_xfm    
    test_jacobian.inputs.jacobianLog = root_dir + log_jacobian
    print test_jacobian.cmdline
#    test_jacobian.run()
