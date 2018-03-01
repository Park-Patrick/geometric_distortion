from nipype.interfaces.base import (
    TraitedSpec,
    CommandLine,
    CommandLineInputSpec,
    File,
    traits
)
import os

class reg_resampleInputSpec(CommandLineInputSpec):
    reference_img = File(argstr='-ref %s', mandatory=True,
        position=0, desc='reference image')

    floating_img = File(argstr='-flo %s',
        position = 1, desc = 'floating image')

    cpp_img = File(argstr='-cpp %s',
        position = 2, desc = 'cpp image')

    aff_img = File(argstr='-aff %s',
        position = 2, desc = 'aff image')

    resampled_img = File(argstr='-res %s',
        position = 3, desc = 'resampled image')


class reg_resampleOutputSpec(TraitedSpec):
    output_resampled_img = File(desc='output resampled image')

class reg_resample(CommandLine):
    _cmd = 'reg_resample'
    input_spec = reg_resampleInputSpec
    output_spec = reg_resampleOutputSpec
    
    def _list_outputs(self):
            outputs = self.output_spec().get()
            
            output_resampled_img = self.inputs.resampled_img
            outputs['output_resampled_img'] = output_resampled_img

            return outputs

if __name__ == '__main__':

    test_resample = reg_resample()
    test_resample.inputs.reference_img = 'refimg'
    test_resample.inputs.floating_img = 'flo_img'
    #test_resample.inputs.cpp_img = 'cpp_img'
    #test_resample.inputs.aff_img = 'aff_img'
    test_resample.inputs.resampled_img = 'res_img'
    print test_resample.cmdline
#    test_jacobian.run()