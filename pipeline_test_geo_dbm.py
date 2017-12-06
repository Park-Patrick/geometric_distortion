##################
# Import modules #
##################

from nipype.interfaces.base import(
    TraitedSpec,
    CommandLineInputSpec,
    CommandLine,
    File,
    traits
)
import nipype.interfaces.ants as ants
import nipype.interfaces.fsl as fsl
import nipype.interfaces.freesurfer as fsurfer

import wrapper_reg_aladin as w_reg_aladin
import wrapper_reg_f3d as w_reg_f3d
import wrapper_reg_jacobian as w_reg_jacobian

import wrapper_reg_transform as w_reg_transform
import wrapper_new_reg_transform as w_new_reg_transform

import wrapper_run_convertNiftiFLOAT32to4D as w_convert4D

import nipype.interfaces.utility as util

#import nipype engine
import nipype.pipeline.engine as pe

###########################
# Set file name variables #
###########################

#root directory
root_dir = '/home/ROBARTS/ppark/Documents/geo_dbm/data2/'
#root_dir = '/data/geo_dbm/data2/'

nipype_dir = 'nipype_data/'

original_dir = 'original/'

data_dir = 'EPI_V053/'

output_dir = '7T/Processed/7Tto3T/'

#input files (not generated within pipeline)
#despot1 = 'Preop/Despot/despot1_tr90_fa18.nii.gz'
#despot1_mask = 'Preop/Processed/BrainMask/BrainMask.nii.gz'
despot1 = 'despot1_tr90_fa18.nii.gz'
despot1_mask = 'BrainMask.nii.gz'

highres = 'MPRAGE_std_N4.nii.gz'
highres_bet = 'MPRAGE_std_N4_bet.nii.gz'
highres_mask = 'MPRAGE_std_N4_bet_BrainMask.nii.gz'

highres_N4 = 'MPRAGE_std_N4_N4.nii.gz'
highres_N4_bet = 'MPRAGE_std_N4_N4_bet.nii.gz'

susc_hires = 'Susc3D_MagEchoAvg.nii.gz'

#output file names
despot1_nuc = 'despot1_tr90_fa18_nuc.nii.gz'
despot1_nuc_masked = 'despot1_tr90_fa18_nuc_masked.nii.gz'

resampled_despot1 = 'despot1_05mm.nii.gz'
resampled_despot1_mask = 'despot1_05mm_BrainMask.nii.gz'

resampled_highres = 'MPRAGE_std_N4_05mm.nii.gz'
resampled_highres_bet = 'MPRAGE_std_N4_bet_05mm.nii.gz'
resampled_highres_mask = 'MPRAGE_std_N4_05mm_BrainMask.nii.gz'

rigid_xfm = '7T_to_3T_rigid_xfm.txt'
rigid_nii = '7T_to_3T_rigid.nii'

affine_xfm = '7T_to_3T_affine_xfm.txt'
affine_nii = '7T_to_3T_affine.nii'

nlin_cpp = '7T_to_3T_nlin_cpp.nii'
nlin_nii = '7T_to_3T_nlin.nii'

log_jacobian = '7T_to_3T_log_jacobian.nii'
mat_jacobian = '7T_to_3T_mat_jacobian.nii'
jacobian = '7T_to_3T_jacobian.nii'

affine_displacement_field_gz = 'affine_displacement_field.nii.gz'
affine_displacement_field_4d_gz = 'affine_displacement_field_4d.nii.gz'

nlin_deformation_field_gz = 'nlin_deformation_field_gz.nii.gz'

nlin_displacement_field_gz = 'nlin_displacement_field_gz.nii.gz'
nlin_displacement_field_4d_gz = 'nlin_displacement_field_4d.nii.gz'
nlin_displacement_field_4d_only = 'nlin_displacement_field_4d_only.nii.gz'
nlin_displacement_field_4d_only_split = 'nlin_displacement_field_4d_only_split'

nlin_displacement_field_4d_only_split1 = 'nlin_displacement_field_4d_only_split0000.nii.gz'
nlin_displacement_field_4d_only_split2 = 'nlin_displacement_field_4d_only_split0001.nii.gz'
nlin_displacement_field_4d_only_split3 = 'nlin_displacement_field_4d_only_split0002.nii.gz'

nlin_displacement_sqr = 'nlin_displacement_sqr.nii.gz'
nlin_displacement = 'nlin_displacement.nii.gz'

nlin_displacement_1mm = 'nlin_displacement_1mm.nii.gz'

resampleVS = (0.5, 0.5, 0.5)

matlab_runtime_lib = '/usr/bin/MATLAB_Runtime/v93'

################
# Create Nodes #
################

# wraps command N4BiasFieldCorrection #

node_N4BiasFieldCorrection = pe.Node(interface=ants.N4BiasFieldCorrection(), name='node_N4BiasFieldCorrection')
node_N4BiasFieldCorrection.inputs.input_image = root_dir + nipype_dir + despot1
node_N4BiasFieldCorrection.inputs.output_image = root_dir + nipype_dir + despot1_nuc
#default save_bias set to false. Output is output_image
#print node_N4BiasFieldCorrection.interface.cmdline

# wraps command fslmaths with -mul flag #

node_fslmath_multiply = pe.Node(interface=fsl.BinaryMaths(), name = 'node_fslmath_multiply')
# uncomment for individual node testing #
#node_fslmath_multiply.inputs.in_file = root_dir + despot1_nuc                                   #check difference with uncomment
#node_fslmath_multiply.inputs.in_file = node_N4BiasFieldCorrection.inputs.output_image
node_fslmath_multiply.inputs.operand_file = root_dir + nipype_dir + despot1_mask
node_fslmath_multiply.inputs.operation = "mul"
node_fslmath_multiply.inputs.out_file = root_dir + nipype_dir + despot1_nuc_masked
#print node_fslmath_multiply.interface.cmdline

# wraps command fslmaths with -bin flag #

node_fslmath_binarize = pe.Node(interface=fsl.UnaryMaths(), name = 'node_fslmath_binarize')
node_fslmath_binarize.inputs.operation = 'bin'
node_fslmath_binarize.inputs.in_file = root_dir + nipype_dir + highres_bet
node_fslmath_binarize.inputs.out_file = root_dir + nipype_dir + highres_mask 
#print node_fslmath_binarize.interface.cmdline

# wraps command mri_convert #
# mRs = mri_Resample

#highres
node_mRs = pe.Node(interface=fsurfer.Resample(), name = 'node_mRs')
node_mRs.inputs.voxel_size = resampleVS
node_mRs.inputs.in_file = root_dir + nipype_dir + highres
node_mRs.inputs.resampled_file = root_dir + nipype_dir + resampled_highres
#print node_mRs.interface.cmdline

#highres_bet
node_mRs_highres_bet = pe.Node(interface=fsurfer.Resample(), name = 'node_mRs_highres_bet')
node_mRs_highres_bet.inputs.voxel_size = resampleVS
node_mRs_highres_bet.inputs.in_file = root_dir + nipype_dir + highres_bet
node_mRs_highres_bet.inputs.resampled_file = root_dir + nipype_dir + resampled_highres_bet
#print node_mRs_highres_bet.interface.cmdline

#highres_mask
node_mRs_highres_mask = pe.Node(interface=fsurfer.Resample(), name = 'node_mRs_highres_mask')
node_mRs_highres_mask.inputs.voxel_size = resampleVS
node_mRs_highres_mask.inputs.in_file = root_dir + nipype_dir + highres_mask
node_mRs_highres_mask.inputs.resampled_file = root_dir + nipype_dir + resampled_highres_mask
#print node_mRs_highres_mask.interface.cmdline

#despot1
node_mRs_despot1 = pe.Node(interface=fsurfer.Resample(), name = 'node_mRs_despot1')
node_mRs_despot1.inputs.voxel_size = resampleVS
#node_mRs_despot1.inputs.in_file = root_dir + nipype_dir + despot1_nuc_masked
node_mRs_despot1.inputs.resampled_file = root_dir + nipype_dir + resampled_despot1
#print node_mRs_despot1.interface.cmdline

#despot1_mask
node_mRs_despot1_mask = pe.Node(interface=fsurfer.Resample(), name = 'node_mRs_despot1_mask')
node_mRs_despot1_mask.inputs.voxel_size = resampleVS
node_mRs_despot1_mask.inputs.in_file = root_dir + nipype_dir + despot1_mask
node_mRs_despot1_mask.inputs.resampled_file = root_dir + nipype_dir + resampled_despot1_mask
#print node_mRs_despot1_mask.interface.cmdline

#niftyReg rigid aladin
node_reg_aladin = pe.Node(interface = w_reg_aladin.reg_aladin(), name = 'node_reg_aladin')
node_reg_aladin.inputs.reference_img = root_dir + nipype_dir + resampled_despot1
node_reg_aladin.inputs.floating_img = root_dir + nipype_dir + resampled_highres_bet
node_reg_aladin.inputs.option_rigOnly = True
node_reg_aladin.inputs.option_affine = root_dir + nipype_dir + rigid_xfm
node_reg_aladin.inputs.resampled_image = root_dir + nipype_dir + rigid_nii
#print node_reg_aladin.interface.cmdline

#niftyReg non rigid aladin
node_reg_aladin2 = pe.Node(interface=w_reg_aladin.reg_aladin(), name = 'node_reg_aladin2')
node_reg_aladin2.inputs.reference_img = root_dir + nipype_dir + resampled_despot1
node_reg_aladin2.inputs.floating_img = root_dir + nipype_dir + resampled_highres_bet
node_reg_aladin2.inputs.input_affine_trans = root_dir + nipype_dir + rigid_xfm
node_reg_aladin2.inputs.option_affine = root_dir + nipype_dir + affine_xfm
node_reg_aladin2.inputs.resampled_image = root_dir + nipype_dir + affine_nii
#print node_reg_aladin2.interface.cmdline

#niftyReg f3d
node_reg_f3d = pe.Node(interface=w_reg_f3d.reg_f3d(), name = 'node_reg_f3d')
node_reg_f3d.inputs.reference_img = root_dir + nipype_dir + resampled_despot1
node_reg_f3d.inputs.floating_img = root_dir + nipype_dir + resampled_highres_bet
node_reg_f3d.inputs.option_affine = root_dir + nipype_dir + affine_xfm
node_reg_f3d.inputs.control_point_grid = root_dir + nipype_dir + nlin_cpp
node_reg_f3d.inputs.resampled_img = root_dir + nipype_dir + nlin_nii
#print node_reg_f3d.interface.cmdline

#niftyReg jacobian Log
node_jacobian_log = pe.Node(interface=w_reg_jacobian.reg_jacobian(), name = 'node_jacobian_log')
node_jacobian_log.inputs.reference_img = root_dir + nipype_dir + resampled_despot1
node_jacobian_log.inputs.control_point_grid = root_dir + nipype_dir + nlin_cpp
node_jacobian_log.inputs.option_affine = root_dir + nipype_dir + affine_xfm
node_jacobian_log.inputs.jacobianLog = root_dir + nipype_dir + log_jacobian
#print node_jacobian_log.interface.cmdline

#niftyReg jacobian matrix
node_jacobian_matrix = pe.Node(interface=w_reg_jacobian.reg_jacobian(), name = 'node_jacobian_matrix')
node_jacobian_matrix.inputs.reference_img = root_dir + nipype_dir + resampled_despot1
node_jacobian_matrix.inputs.control_point_grid = root_dir + nipype_dir + nlin_cpp
node_jacobian_matrix.inputs.option_affine = root_dir + nipype_dir + affine_xfm    
node_jacobian_matrix.inputs.jacobianMatrix = root_dir + nipype_dir + mat_jacobian
#print node_jacobian_matrix.interface.cmdline

#niftyReg jacobian determinant
node_jacobian_dtm = pe.Node(interface=w_reg_jacobian.reg_jacobian(), name = 'node_jacobian_dtm')
node_jacobian_dtm.inputs.reference_img = root_dir + nipype_dir + resampled_despot1
node_jacobian_dtm.inputs.control_point_grid = root_dir + nipype_dir +  nlin_cpp
node_jacobian_dtm.inputs.option_affine = root_dir + nipype_dir + affine_xfm    
node_jacobian_dtm.inputs.jacobianDeterminant = root_dir + nipype_dir + jacobian
#print node_jacobian_dtm.interface.cmdline

#niftyReg new transform
node_new_reg_transform = pe.Node(interface=w_new_reg_transform.new_reg_transform(), name = 'node_new_reg_transform')
node_new_reg_transform.inputs.reference_img = root_dir + nipype_dir + resampled_despot1
node_new_reg_transform.inputs.disp = root_dir + nipype_dir + affine_xfm
node_new_reg_transform.inputs.disp2 = root_dir + nipype_dir + affine_displacement_field_gz

#niftyReg cpp transform
node_cpp_reg_transform = pe.Node(interface=w_reg_transform.reg_transform(), name = 'node_cpp_reg_transform')
node_cpp_reg_transform.inputs.reference_img = root_dir + nipype_dir + resampled_despot1
node_cpp_reg_transform.inputs.cpp2def = root_dir + nipype_dir + nlin_cpp
node_cpp_reg_transform.inputs.cpp2def_2 = root_dir + nipype_dir + nlin_deformation_field_gz

#niftyReg def transform
node_def_reg_transform = pe.Node(interface=w_reg_transform.reg_transform(), name = 'node_def_reg_transform')
node_def_reg_transform.inputs.reference_img = root_dir + nipype_dir + resampled_despot1
node_def_reg_transform.inputs.def2disp = root_dir + nipype_dir + nlin_deformation_field_gz
node_def_reg_transform.inputs.def2disp_2 = root_dir + nipype_dir + nlin_displacement_field_gz

#affine 4d convert
node_aff_4d_convert = pe.Node(interface=w_convert4D.run_convertNiftiFLOAT32to4D(), name = 'node_aff_4d_convert')
node_aff_4d_convert.inputs.runtime_lib = matlab_runtime_lib
node_aff_4d_convert.inputs.in_nii = root_dir + nipype_dir + affine_displacement_field_gz
node_aff_4d_convert.inputs.out_4d_nii = root_dir + nipype_dir + affine_displacement_field_4d_gz

#nlin 4d convert
node_nlin_4d_convert = pe.Node(interface=w_convert4D.run_convertNiftiFLOAT32to4D(), name = 'node_nlin_4d_convert')
node_nlin_4d_convert.inputs.runtime_lib = matlab_runtime_lib
node_nlin_4d_convert.inputs.in_nii = root_dir + nipype_dir + nlin_displacement_field_gz
node_nlin_4d_convert.inputs.out_4d_nii = root_dir + nipype_dir + nlin_displacement_field_4d_gz

#fsl subtract
node_fsl_subtract = pe.Node(interface=fsl.BinaryMaths(), name = 'node_fsl_subtract')
#node_fsl_subtract.inputs.in_file = root_dir + nipype_dir + nlin_displacement_field_4d_gz
#node_fsl_subtract.inputs.operand_file = root_dir + nipype_dir + affine_displacement_field_4d_gz
node_fsl_subtract.inputs.operation = 'sub'
node_fsl_subtract.inputs.out_file = root_dir + nipype_dir + nlin_displacement_field_4d_only

#fsl split
node_fsl_split = pe.Node(interface=fsl.Split(), name = 'node_fsl_split')
node_fsl_split.inputs.dimension = 't'
#node_fsl_split.inputs.in_file = root_dir + nipype_dir + nlin_displacement_field_4d_only
node_fsl_split.inputs.out_base_name = root_dir + nipype_dir + nlin_displacement_field_4d_only_split

#select node1
node_select1 = pe.Node(interface=util.Select(), name = 'node_select1')
node_select1.inputs.index=[0]
node_select1.inputs.inlist=[root_dir + nipype_dir + nlin_displacement_field_4d_only_split1,
							root_dir + nipype_dir + nlin_displacement_field_4d_only_split2,
							root_dir + nipype_dir + nlin_displacement_field_4d_only_split3]

#select node2
node_select2 = pe.Node(interface=util.Select(), name = 'node_select2')
node_select2.inputs.index=[1]
node_select2.inputs.inlist=[root_dir + nipype_dir + nlin_displacement_field_4d_only_split1,
							root_dir + nipype_dir + nlin_displacement_field_4d_only_split2,
							root_dir + nipype_dir + nlin_displacement_field_4d_only_split3]

#select node3
node_select3 = pe.Node(interface=util.Select(), name = 'node_select3')
node_select3.inputs.index=[2]
node_select3.inputs.inlist=[root_dir + nipype_dir + nlin_displacement_field_4d_only_split1,
							root_dir + nipype_dir + nlin_displacement_field_4d_only_split2,
							root_dir + nipype_dir + nlin_displacement_field_4d_only_split3]

#fsl square 1
node_fsl_sqr1 = pe.Node(interface=fsl.UnaryMaths(), name = 'node_fsl_sqr1')
node_fsl_sqr1.inputs.operation = 'sqr'
#node_fsl_sqr1.inputs.in_file = root_dir + nipype_dir + nlin_displacement_field_4d_only_split1
node_fsl_sqr1.inputs.out_file = root_dir + nipype_dir + nlin_displacement_field_4d_only_split1
#print node_fsl_sqr1.interface.cmdline

#fsl square 2
node_fsl_sqr2 = pe.Node(interface=fsl.UnaryMaths(), name = 'node_fsl_sqr2')
node_fsl_sqr2.inputs.operation = 'sqr'
#node_fsl_sqr2.inputs.in_file = root_dir + nipype_dir + nlin_displacement_field_4d_only_split2
node_fsl_sqr2.inputs.out_file = root_dir + nipype_dir + nlin_displacement_field_4d_only_split2
#print node_fsl_sqr2.interface.cmdline

#fsl square 3
node_fsl_sqr3 = pe.Node(interface=fsl.UnaryMaths(), name = 'node_fsl_sqr3')
node_fsl_sqr3.inputs.operation = 'sqr'
#node_fsl_sqr3.inputs.in_file = root_dir + nipype_dir + nlin_displacement_field_4d_only_split3
node_fsl_sqr3.inputs.out_file = root_dir + nipype_dir + nlin_displacement_field_4d_only_split3
#print node_fsl_sqr3.interface.cmdline

#merge node
node_merge = pe.Node(interface=util.Merge(2), name = 'node_merge')
node_merge.inputs.in1 = root_dir + nipype_dir + nlin_displacement_field_4d_only_split2
node_merge.inputs.in2 = root_dir + nipype_dir + nlin_displacement_field_4d_only_split3

#fsl add all
node_fsl_add_all = pe.Node(interface=fsl.MultiImageMaths(), name = 'node_fsl_add_all')
#node_fsl_add_all.inputs.in_file = root_dir + nipype_dir + nlin_displacement_field_4d_only_split1
node_fsl_add_all.inputs.op_string = '-add %s -add %s'
#node_fsl_add_all.inputs.operand_files = [root_dir + nipype_dir + nlin_displacement_field_4d_only_split2, root_dir + nipype_dir + nlin_displacement_field_4d_only_split3]
node_fsl_add_all.inputs.out_file = root_dir + nipype_dir + nlin_displacement_sqr

#fsl square root
node_fsl_sqrt = pe.Node(interface=fsl.UnaryMaths(), name = 'node_fsl_sqrt')
#node_fsl_sqrt.inputs.in_file = root_dir + nipype_dir + nlin_displacement_sqr
node_fsl_sqrt.inputs.operation = 'sqrt'
node_fsl_sqrt.inputs.out_file = root_dir + nipype_dir + nlin_displacement

#reg resample

####################
# Create workflows #
####################
#TODO: break pipeline down into sub-workflows

test_geo_workflow = pe.Workflow(name='test_geo_workflow')
test_geo_workflow.base_dir = root_dir + 'base/'

#############
# Add nodes #
#############
#TODO: can be simplified to one line per sub-workflow


test_geo_workflow.add_nodes([node_fslmath_binarize])
test_geo_workflow.add_nodes([node_fslmath_multiply])
test_geo_workflow.add_nodes([node_jacobian_dtm])
test_geo_workflow.add_nodes([node_jacobian_log])
test_geo_workflow.add_nodes([node_jacobian_matrix])
test_geo_workflow.add_nodes([node_mRs])
test_geo_workflow.add_nodes([node_mRs_despot1])
test_geo_workflow.add_nodes([node_mRs_despot1_mask])
test_geo_workflow.add_nodes([node_mRs_highres_bet])
test_geo_workflow.add_nodes([node_mRs_highres_mask])
test_geo_workflow.add_nodes([node_N4BiasFieldCorrection])
test_geo_workflow.add_nodes([node_reg_aladin])
test_geo_workflow.add_nodes([node_reg_aladin2])
test_geo_workflow.add_nodes([node_reg_f3d])
test_geo_workflow.add_nodes([node_new_reg_transform])
test_geo_workflow.add_nodes([node_cpp_reg_transform])
test_geo_workflow.add_nodes([node_def_reg_transform])
test_geo_workflow.add_nodes([node_aff_4d_convert])
test_geo_workflow.add_nodes([node_nlin_4d_convert])
test_geo_workflow.add_nodes([node_fsl_subtract])
test_geo_workflow.add_nodes([node_fsl_split])
test_geo_workflow.add_nodes([node_select1])
test_geo_workflow.add_nodes([node_select2])
test_geo_workflow.add_nodes([node_select3])
test_geo_workflow.add_nodes([node_fsl_sqr1])
test_geo_workflow.add_nodes([node_fsl_sqr2])
test_geo_workflow.add_nodes([node_fsl_sqr3])
test_geo_workflow.add_nodes([node_merge])
test_geo_workflow.add_nodes([node_fsl_add_all])
test_geo_workflow.add_nodes([node_fsl_sqrt])

#################
# Connect nodes #
#################
#TODO: can be simplified to one line per sub-workflow

test_geo_workflow.connect(node_N4BiasFieldCorrection, 'output_image', node_fslmath_multiply, 'in_file')
test_geo_workflow.connect(node_fslmath_multiply, 'out_file', node_mRs_despot1, 'in_file')
test_geo_workflow.connect(node_mRs_highres_bet, 'resampled_file', node_reg_aladin, 'floating_img')
test_geo_workflow.connect(node_mRs_highres_bet, 'resampled_file', node_reg_aladin2, 'floating_img')
test_geo_workflow.connect(node_reg_aladin, 'output_affine_transformation', node_reg_aladin2, 'input_affine_trans')
#reference images
test_geo_workflow.connect(node_mRs_despot1, 'resampled_file', node_reg_aladin, 'reference_img')
test_geo_workflow.connect(node_mRs_despot1, 'resampled_file', node_reg_aladin2, 'reference_img')
test_geo_workflow.connect(node_mRs_despot1, 'resampled_file', node_reg_f3d, 'reference_img')
test_geo_workflow.connect(node_mRs_despot1, 'resampled_file', node_jacobian_dtm, 'reference_img')
test_geo_workflow.connect(node_mRs_despot1, 'resampled_file', node_jacobian_log, 'reference_img')
test_geo_workflow.connect(node_mRs_despot1, 'resampled_file', node_jacobian_matrix, 'reference_img')

test_geo_workflow.connect(node_mRs_highres_bet, 'resampled_file', node_reg_f3d, 'floating_img')
test_geo_workflow.connect(node_reg_aladin2, 'output_affine_transformation', node_reg_f3d, 'option_affine')

test_geo_workflow.connect(node_reg_f3d, 'outfile_cpp', node_jacobian_log, 'control_point_grid')
test_geo_workflow.connect(node_reg_f3d, 'outfile_cpp', node_jacobian_matrix, 'control_point_grid')
test_geo_workflow.connect(node_reg_f3d, 'outfile_cpp', node_jacobian_dtm, 'control_point_grid')

test_geo_workflow.connect(node_reg_aladin2, 'output_affine_transformation', node_jacobian_log, 'option_affine')
test_geo_workflow.connect(node_reg_aladin2, 'output_affine_transformation', node_jacobian_matrix, 'option_affine')
test_geo_workflow.connect(node_reg_aladin2, 'output_affine_transformation', node_jacobian_dtm, 'option_affine')

test_geo_workflow.connect(node_mRs_despot1, 'resampled_file', node_new_reg_transform, 'reference_img')
test_geo_workflow.connect(node_reg_aladin2, 'output_affine_transformation', node_new_reg_transform, 'disp')

test_geo_workflow.connect(node_mRs_despot1, 'resampled_file', node_cpp_reg_transform, 'reference_img')
test_geo_workflow.connect(node_reg_f3d, 'outfile_cpp', node_cpp_reg_transform, 'cpp2def')

test_geo_workflow.connect(node_mRs_despot1, 'resampled_file', node_def_reg_transform, 'reference_img')
test_geo_workflow.connect(node_cpp_reg_transform, 'output_nlin_deformation', node_def_reg_transform, 'def2disp')

test_geo_workflow.connect(node_new_reg_transform, 'output_disp2', node_aff_4d_convert, 'in_nii')
test_geo_workflow.connect(node_def_reg_transform, 'output_nlin_displacement', node_nlin_4d_convert, 'in_nii')

test_geo_workflow.connect(node_nlin_4d_convert, 'output_4d_nii', node_fsl_subtract, 'in_file')
test_geo_workflow.connect(node_aff_4d_convert, 'output_4d_nii', node_fsl_subtract, 'operand_file')

test_geo_workflow.connect(node_fsl_subtract, 'out_file', node_fsl_split, 'in_file')

test_geo_workflow.connect(node_fsl_split, 'out_files', node_select1, 'inlist')
test_geo_workflow.connect(node_select1, 'out', node_fsl_sqr1, 'in_file')

test_geo_workflow.connect(node_fsl_split, 'out_files', node_select2, 'inlist')
test_geo_workflow.connect(node_select2, 'out', node_fsl_sqr2, 'in_file')

test_geo_workflow.connect(node_fsl_split, 'out_files', node_select3, 'inlist')
test_geo_workflow.connect(node_select3, 'out', node_fsl_sqr3, 'in_file')

test_geo_workflow.connect(node_fsl_sqr2, 'out_file', node_merge, 'in1')
test_geo_workflow.connect(node_fsl_sqr3, 'out_file', node_merge, 'in2')

test_geo_workflow.connect(node_fsl_sqr1, 'out_file', node_fsl_add_all, 'in_file') 
test_geo_workflow.connect(node_merge, 'out', node_fsl_add_all, 'operand_files')

test_geo_workflow.connect(node_fsl_add_all, 'out_file', node_fsl_sqrt, 'in_file')
##############
# Draw graph #
##############

#simple
test_geo_workflow.write_graph()

#detailed
test_geo_workflow.write_graph('test_geo_workflow_graph.dot', graph2use = 'exec')

#######
# Run #
#######


test_geo_workflow.run()




