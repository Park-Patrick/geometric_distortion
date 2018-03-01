from nipype.interfaces.base import (
    TraitedSpec,
    CommandLine,
    CommandLineInputSpec,
    File,
    traits
)
import os

class reg_transformInputSpec(CommandLineInputSpec):
    reference_img = File(argstr='-ref %s', mandatory=True,
                         position=0, desc='reference image')

    cpp2def = File(argstr='-cpp2def %s',
                              position = 1, desc = 'cpp2def')

    cpp2def_2 = File(argstr='%s',
                            position = 2, desc = 'cpp2def_2')

    def2disp = File(argstr='-def2disp %s',
                         position = 1, desc = 'def2disp')

    def2disp_2 = File(argstr='%s',
                            position = 2, desc = 'def2disp_2')
                         

class reg_transformOutputSpec(TraitedSpec):
    output_nlin_deformation = File(desc='output nlin deformation')
    output_nlin_displacement = File(desc='output nlin displacement')

class reg_transform(CommandLine):
    _cmd = 'reg_transform'
    input_spec = reg_transformInputSpec
    output_spec = reg_transformOutputSpec
    
    def _list_outputs(self):
            outputs = self.output_spec().get()
            
            output_nlin_deformation = self.inputs.cpp2def_2
            outputs['output_nlin_deformation'] = output_nlin_deformation

            output_nlin_displacement = self.inputs.def2disp_2
            outputs['output_nlin_displacement'] = output_nlin_displacement

            return outputs

if __name__ == '__main__':

    test_transform = reg_transform()
    test_transform.inputs.reference_img = 'refimg'
    #test_transform.inputs.cpp2def = 'cpp'
    #test_transform.inputs.cpp2def_2 = 'def'
    test_transform.inputs.def2disp = 'def'
    test_transform.inputs.def2disp_2 = 'disp'
    print test_transform.cmdline
    test_transform.run()
