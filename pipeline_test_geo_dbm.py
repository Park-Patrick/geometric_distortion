

## Pipeline wrapper code

import os
import sys

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
import wrapper_reg_resample as w_reg_resample

import wrapper_run_convertNiftiFLOAT32to4D as w_convert4D

import nipype.interfaces.utility as util

#import nipype engine
import nipype.pipeline.engine as pe

subj_id = sys.argv[1]

DIS2_or_3 = sys.argv[2]

#root directory
#root_dir = '/home/ROBARTS/ppark/Documents/geo_dbm/data3/'
root_dir = '/data/OHBM_dbm/'

nipype_dir = 'sub-' + subj_id + "/"

# original_dir = 'original/'

# data_dir = 'EPI_V053/'

# output_dir = '7T/Processed/7Tto3T/'

workflow_dir = nipype_dir

#######################
##### OHBM inputs #####
#######################

#resampled to --> resampled_despot1
despot1_nuc_masked = 'sub-' + subj_id + '_acq-MP2RAGE_rec-' + DIS2_or_3 + '_T1w.nii.gz'

#resampled to --> resampled_highres_bet
highres_N4_bet = 'sub-' + subj_id + '_acq-MP2RAGE_T1w.nii.gz'


#input files (not generated within pipeline)
#despot1 = 'Preop/Despot/despot1_tr90_fa18.nii.gz'
#despot1_mask = 'Preop/Processed/BrainMask/BrainMask.nii.gz'
#despot1 = 'despot1_tr90_fa18.nii.gz'
#despot1_mask = 'BrainMask.nii.gz'

#highres = 'MPRAGE_std_N4.nii.gz'
#highres_bet = 'MPRAGE_std_N4_bet.nii.gz'
#highres_mask = 'MPRAGE_std_N4_bet_BrainMask.nii.gz'

#susc_hires = 'Susc3D_MagEchoAvg.nii.gz'

#output file names
#despot1_nuc = 'despot1_tr90_fa18_nuc.nii.gz'
#despot1_nuc_masked = 'despot1_tr90_fa18_nuc_masked.nii.gz'

#highres_N4 = 'MPRAGE_std_N4_N4.nii.gz'
#highres_N4_bet = 'MPRAGE_std_N4_N4_bet.nii.gz'

resampled_despot1 = 'despot1_05mm.nii.gz'
#resampled_despot1_mask = 'despot1_05mm_BrainMask.nii.gz'

#resampled_highres = 'MPRAGE_std_N4_05mm.nii.gz'
resampled_highres_bet = 'MPRAGE_std_N4_bet_05mm.nii.gz'
#resampled_highres_mask = 'MPRAGE_std_N4_05mm_BrainMask.nii.gz'


rigid_xfm_7T_to_3T = '7T_to_3T_rigid_xfm.txt'
rigid_nii_7T_to_3T = '7T_to_3T_rigid.nii'

affine_xfm_7T_to_3T = '7T_to_3T_affine_xfm.txt'
affine_nii_7T_to_3T = '7T_to_3T_affine.nii'

nlin_cpp_7T_to_3T = '7T_to_3T_nlin_cpp.nii'
nlin_nii_7T_to_3T = '7T_to_3T_nlin.nii'

log_jacobian_7T_to_3T = '7T_to_3T_log_jacobian.nii'
mat_jacobian_7T_to_3T = '7T_to_3T_mat_jacobian.nii'
jacobian_7T_to_3T = '7T_to_3T_jacobian.nii'


rigid_xfm_3T_to_MNI = '3T_to_MNI_rigid_xfm.txt'
rigid_nii_3T_to_MNI= '3T_to_MNI_rigid.nii'

affine_xfm_3T_to_MNI = '3T_to_MNI_affine_xfm.txt'
affine_nii_3T_to_MNI = '3T_to_MNI_affine.nii'

nlin_cpp_3T_to_MNI = '3T_to_MNI_nlin_cpp.nii'
nlin_nii_3T_to_MNI = '3T_to_MNI_nlin.nii'

log_jacobian_3T_to_MNI= '3T_to_MNI_log_jacobian.nii'
mat_jacobian_3T_to_MNI = '3T_to_MNI_mat_jacobian.nii'
jacobian_3T_to_MNI = '3T_to_MNI_jacobian.nii'


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

x_nlin_displacement = 'nlin_displacement_x.nii.gz'
y_nlin_displacement = 'nlin_displacement_y.nii.gz'
z_nlin_displacement = 'nlin_displacement_z.nii.gz'

resampled_7T_to_MNI = '7T_to_MNI_nlin.nii'

resampled_atlas = root_dir + 'atlas/t1.05mm.nii.gz'
resampled_atlas_1mm= root_dir + 'atlas/t1.nii.gz'
resampled_atlas_masked = root_dir + 'atlas/t1.brain.inorm.05mm.nii.gz'

subjMNI_nlin_displacement_x_1mm = 'subj_MNI_displacement_nlin_x_1mm.nii.gz'
subjMNI_nlin_displacement_y_1mm = 'subj_MNI_displacement_nlin_y_1mm.nii.gz'
subjMNI_nlin_displacement_z_1mm = 'subj_MNI_displacement_nlin_z_1mm.nii.gz'

subjMNI_nlin_displacement = 'subj_MNI_displacement_nlin.nii.gz'
subjMNI_nlin_displacement_1mm = 'subj_MNI_displacement_nlin_1mm.nii.gz'

resampleVS = (0.5, 0.5, 0.5)

matlab_runtime_lib = '/usr/bin/MATLAB_Runtime/v93'


# Node 24
node_mRs_highres_bet = pe.Node(interface=fsurfer.Resample(), name = 'node_mRs_highres_bet')
node_mRs_highres_bet.inputs.voxel_size = resampleVS
node_mRs_highres_bet.inputs.in_file = root_dir + nipype_dir + highres_N4_bet
node_mRs_highres_bet.inputs.resampled_file = root_dir + nipype_dir + resampled_highres_bet

# Node 25
node_mRs_despot1 = pe.Node(interface=fsurfer.Resample(), name = 'node_mRs_despot1')
node_mRs_despot1.inputs.voxel_size = resampleVS
node_mRs_despot1.inputs.in_file = root_dir + nipype_dir + despot1_nuc_masked
node_mRs_despot1.inputs.resampled_file = root_dir + nipype_dir + resampled_despot1

# Node 1
node_reg_aladin = pe.Node(interface = w_reg_aladin.reg_aladin(), name = 'node_reg_aladin')
node_reg_aladin.inputs.reference_img = root_dir + nipype_dir + resampled_despot1
node_reg_aladin.inputs.floating_img = root_dir + nipype_dir + resampled_highres_bet
node_reg_aladin.inputs.option_rigOnly = True
node_reg_aladin.inputs.option_affine = root_dir + nipype_dir + rigid_xfm_7T_to_3T
node_reg_aladin.inputs.resampled_image = root_dir + nipype_dir + rigid_nii_7T_to_3T

# Node 2
node_reg_aladin2 = pe.Node(interface=w_reg_aladin.reg_aladin(), name = 'node_reg_aladin2')
node_reg_aladin2.inputs.reference_img = root_dir + nipype_dir + resampled_despot1
node_reg_aladin2.inputs.floating_img = root_dir + nipype_dir + resampled_highres_bet
node_reg_aladin2.inputs.input_affine_trans = root_dir + nipype_dir + rigid_xfm_7T_to_3T
node_reg_aladin2.inputs.option_affine = root_dir + nipype_dir + affine_xfm_7T_to_3T
node_reg_aladin2.inputs.resampled_image = root_dir + nipype_dir + affine_nii_7T_to_3T

# Node 3
node_reg_f3d = pe.Node(interface=w_reg_f3d.reg_f3d(), name = 'node_reg_f3d')
node_reg_f3d.inputs.reference_img = root_dir + nipype_dir + resampled_despot1
node_reg_f3d.inputs.floating_img = root_dir + nipype_dir + resampled_highres_bet
node_reg_f3d.inputs.option_affine = root_dir + nipype_dir + affine_xfm_7T_to_3T
node_reg_f3d.inputs.control_point_grid = root_dir + nipype_dir + nlin_cpp_7T_to_3T
node_reg_f3d.inputs.resampled_img = root_dir + nipype_dir + nlin_nii_7T_to_3T

# Node 4
node_new_reg_transform = pe.Node(interface=w_new_reg_transform.new_reg_transform(), name = 'node_new_reg_transform')
node_new_reg_transform.inputs.reference_img = root_dir + nipype_dir + resampled_despot1
node_new_reg_transform.inputs.disp = root_dir + nipype_dir + affine_xfm_7T_to_3T
node_new_reg_transform.inputs.disp2 = root_dir + nipype_dir + affine_displacement_field_gz

# Node 5
node_cpp_reg_transform = pe.Node(interface=w_reg_transform.reg_transform(), name = 'node_cpp_reg_transform')
node_cpp_reg_transform.inputs.reference_img = root_dir + nipype_dir + resampled_despot1
node_cpp_reg_transform.inputs.cpp2def = root_dir + nipype_dir + nlin_cpp_7T_to_3T
node_cpp_reg_transform.inputs.cpp2def_2 = root_dir + nipype_dir + nlin_deformation_field_gz

# Node 6
node_aff_4d_convert = pe.Node(interface=w_convert4D.run_convertNiftiFLOAT32to4D(), name = 'node_aff_4d_convert')
node_aff_4d_convert.inputs.runtime_lib = matlab_runtime_lib
node_aff_4d_convert.inputs.in_nii = root_dir + nipype_dir + affine_displacement_field_gz
node_aff_4d_convert.inputs.out_4d_nii = root_dir + nipype_dir + affine_displacement_field_4d_gz

# Node 7
node_def_reg_transform = pe.Node(interface=w_reg_transform.reg_transform(), name = 'node_def_reg_transform')
node_def_reg_transform.inputs.reference_img = root_dir + nipype_dir + resampled_despot1
node_def_reg_transform.inputs.def2disp = root_dir + nipype_dir + nlin_deformation_field_gz
node_def_reg_transform.inputs.def2disp_2 = root_dir + nipype_dir + nlin_displacement_field_gz

# Node 8
node_nlin_4d_convert = pe.Node(interface=w_convert4D.run_convertNiftiFLOAT32to4D(), name = 'node_nlin_4d_convert')
node_nlin_4d_convert.inputs.runtime_lib = matlab_runtime_lib
node_nlin_4d_convert.inputs.in_nii = root_dir + nipype_dir + nlin_displacement_field_gz
node_nlin_4d_convert.inputs.out_4d_nii = root_dir + nipype_dir + nlin_displacement_field_4d_gz

# Node 9
node_fsl_subtract = pe.Node(interface=fsl.BinaryMaths(), name = 'node_fsl_subtract')
#node_fsl_subtract.inputs.in_file = root_dir + nipype_dir + nlin_displacement_field_4d_gz
#node_fsl_subtract.inputs.operand_file = root_dir + nipype_dir + affine_displacement_field_4d_gz
node_fsl_subtract.inputs.operation = 'sub'
node_fsl_subtract.inputs.out_file = root_dir + nipype_dir + nlin_displacement_field_4d_only

# Node 10
node_fsl_split = pe.Node(interface=fsl.Split(), name = 'node_fsl_split')
node_fsl_split.inputs.dimension = 't'
#node_fsl_split.inputs.in_file = root_dir + nipype_dir + nlin_displacement_field_4d_only
node_fsl_split.inputs.out_base_name = root_dir + nipype_dir + nlin_displacement_field_4d_only_split

# Node 11
node_select1 = pe.Node(interface=util.Select(), name = 'node_select1')
node_select1.inputs.index=[0]
node_select1.inputs.inlist=[root_dir + nipype_dir + nlin_displacement_field_4d_only_split1,
							root_dir + nipype_dir + nlin_displacement_field_4d_only_split2,
							root_dir + nipype_dir + nlin_displacement_field_4d_only_split3]

# Node 12
node_select2 = pe.Node(interface=util.Select(), name = 'node_select2')
node_select2.inputs.index=[1]
node_select2.inputs.inlist=[root_dir + nipype_dir + nlin_displacement_field_4d_only_split1,
							root_dir + nipype_dir + nlin_displacement_field_4d_only_split2,
							root_dir + nipype_dir + nlin_displacement_field_4d_only_split3]

# Node 13
node_select3 = pe.Node(interface=util.Select(), name = 'node_select3')
node_select3.inputs.index=[2]
node_select3.inputs.inlist=[root_dir + nipype_dir + nlin_displacement_field_4d_only_split1,
							root_dir + nipype_dir + nlin_displacement_field_4d_only_split2,
							root_dir + nipype_dir + nlin_displacement_field_4d_only_split3]

# Node 14
node_fsl_sqr1 = pe.Node(interface=fsl.UnaryMaths(), name = 'node_fsl_sqr1')
node_fsl_sqr1.inputs.operation = 'sqr'
#node_fsl_sqr1.inputs.in_file = root_dir + nipype_dir + nlin_displacement_field_4d_only_split1
node_fsl_sqr1.inputs.out_file = root_dir + nipype_dir + nlin_displacement_field_4d_only_split1

# Node 15
node_fsl_sqr2 = pe.Node(interface=fsl.UnaryMaths(), name = 'node_fsl_sqr2')
node_fsl_sqr2.inputs.operation = 'sqr'
#node_fsl_sqr2.inputs.in_file = root_dir + nipype_dir + nlin_displacement_field_4d_only_split2
node_fsl_sqr2.inputs.out_file = root_dir + nipype_dir + nlin_displacement_field_4d_only_split2

# Node 16
node_fsl_sqr3 = pe.Node(interface=fsl.UnaryMaths(), name = 'node_fsl_sqr3')
node_fsl_sqr3.inputs.operation = 'sqr'
#node_fsl_sqr3.inputs.in_file = root_dir + nipype_dir + nlin_displacement_field_4d_only_split3
node_fsl_sqr3.inputs.out_file = root_dir + nipype_dir + nlin_displacement_field_4d_only_split3

# Node 17
node_merge = pe.Node(interface=util.Merge(2), name = 'node_merge')
node_merge.inputs.in1 = root_dir + nipype_dir + nlin_displacement_field_4d_only_split2
node_merge.inputs.in2 = root_dir + nipype_dir + nlin_displacement_field_4d_only_split3

# Node 18
node_fsl_add_all = pe.Node(interface=fsl.MultiImageMaths(), name = 'node_fsl_add_all')
#node_fsl_add_all.inputs.in_file = root_dir + nipype_dir + nlin_displacement_field_4d_only_split1
node_fsl_add_all.inputs.op_string = '-add %s -add %s'
#node_fsl_add_all.inputs.operand_files = [root_dir + nipype_dir + nlin_displacement_field_4d_only_split2, root_dir + nipype_dir + nlin_displacement_field_4d_only_split3]
node_fsl_add_all.inputs.out_file = root_dir + nipype_dir + nlin_displacement_sqr

# Node 19
node_fsl_sqrt = pe.Node(interface=fsl.UnaryMaths(), name = 'node_fsl_sqrt')
#node_fsl_sqrt.inputs.in_file = root_dir + nipype_dir + nlin_displacement_sqr
node_fsl_sqrt.inputs.operation = 'sqrt'
node_fsl_sqrt.inputs.out_file = root_dir + nipype_dir + nlin_displacement

# Node 20
node_3T_MNI_reg_aladin1 = pe.Node(interface = w_reg_aladin.reg_aladin(), name = 'node_3T_MNI_reg_aladin1')
node_3T_MNI_reg_aladin1.inputs.reference_img = resampled_atlas_masked
node_3T_MNI_reg_aladin1.inputs.floating_img = root_dir + nipype_dir + resampled_despot1
node_3T_MNI_reg_aladin1.inputs.option_rigOnly = True
node_3T_MNI_reg_aladin1.inputs.option_affine = root_dir + nipype_dir + rigid_xfm_3T_to_MNI
node_3T_MNI_reg_aladin1.inputs.resampled_image = root_dir + nipype_dir + rigid_nii_3T_to_MNI

# Node 21
node_3T_MNI_reg_aladin2 = pe.Node(interface = w_reg_aladin.reg_aladin(), name = 'node_3T_MNI_reg_aladin2')
node_3T_MNI_reg_aladin2.inputs.reference_img = resampled_atlas_masked
node_3T_MNI_reg_aladin2.inputs.floating_img = root_dir + nipype_dir + resampled_despot1
node_3T_MNI_reg_aladin2.inputs.input_affine_trans = root_dir + nipype_dir + rigid_xfm_3T_to_MNI
node_3T_MNI_reg_aladin2.inputs.option_affine = root_dir + nipype_dir + affine_xfm_3T_to_MNI
node_3T_MNI_reg_aladin2.inputs.resampled_image = root_dir + nipype_dir + affine_nii_3T_to_MNI

# Node 22
node_3T_MNI_reg_f3d = pe.Node(interface=w_reg_f3d.reg_f3d(), name = 'node_3T_MNI_reg_f3d')
node_3T_MNI_reg_f3d.inputs.reference_img = resampled_atlas_masked
node_3T_MNI_reg_f3d.inputs.floating_img = root_dir + nipype_dir + resampled_despot1
node_3T_MNI_reg_f3d.inputs.option_affine = root_dir + nipype_dir + affine_xfm_3T_to_MNI
node_3T_MNI_reg_f3d.inputs.control_point_grid = root_dir + nipype_dir + nlin_cpp_3T_to_MNI
node_3T_MNI_reg_f3d.inputs.resampled_img = root_dir + nipype_dir + nlin_nii_3T_to_MNI

# Node 23-1
node_subjMNI_resample = pe.Node(interface=w_reg_resample.reg_resample(), name = 'node_subjMNI_resample')
node_subjMNI_resample.inputs.reference_img = resampled_atlas
node_subjMNI_resample.inputs.floating_img = root_dir + nipype_dir + nlin_displacement
node_subjMNI_resample.inputs.cpp_img = root_dir + nipype_dir + nlin_cpp_3T_to_MNI
node_subjMNI_resample.inputs.resampled_img = root_dir + nipype_dir + subjMNI_nlin_displacement

# Node 23-2
node_1mm_subjMNI_resample = pe.Node(interface=w_reg_resample.reg_resample(), name = 'node_1mm_subjMNI_resample')
node_1mm_subjMNI_resample.inputs.reference_img = resampled_atlas_1mm
node_1mm_subjMNI_resample.inputs.floating_img = root_dir + nipype_dir + nlin_displacement
node_1mm_subjMNI_resample.inputs.cpp_img = root_dir + nipype_dir + nlin_cpp_3T_to_MNI
node_1mm_subjMNI_resample.inputs.resampled_img = root_dir + nipype_dir + subjMNI_nlin_displacement_1mm

####################
# Create workflows #
####################
#TODO: break pipeline down into sub-workflows

OHBM_workflow = pe.Workflow(name='OHBM_workflow')
OHBM_workflow.base_dir = root_dir + workflow_dir

# Node 1
OHBM_workflow.add_nodes([node_reg_aladin])

# Node 2
OHBM_workflow.add_nodes([node_reg_aladin2])

# Node 3
OHBM_workflow.add_nodes([node_reg_f3d])

# Node 4
OHBM_workflow.add_nodes([node_new_reg_transform])

# Node 5
OHBM_workflow.add_nodes([node_cpp_reg_transform])

# Node 6
OHBM_workflow.add_nodes([node_aff_4d_convert])

# Node 7
OHBM_workflow.add_nodes([node_def_reg_transform])

# Node 8
OHBM_workflow.add_nodes([node_nlin_4d_convert])

# Node 9
OHBM_workflow.add_nodes([node_fsl_subtract])

# Node 10
OHBM_workflow.add_nodes([node_fsl_split])

# Node 11
OHBM_workflow.add_nodes([node_select1])

# Node 12
OHBM_workflow.add_nodes([node_select2])

# Node 13
OHBM_workflow.add_nodes([node_select3])

# Node 14
OHBM_workflow.add_nodes([node_fsl_sqr1])

# Node 15
OHBM_workflow.add_nodes([node_fsl_sqr2])

# Node 16
OHBM_workflow.add_nodes([node_fsl_sqr3])

# Node 17
OHBM_workflow.add_nodes([node_merge])

# Node 18
OHBM_workflow.add_nodes([node_fsl_add_all])

# Node 19
OHBM_workflow.add_nodes([node_fsl_sqrt])

# Node 20
OHBM_workflow.add_nodes([node_3T_MNI_reg_aladin1])

# Node 21
OHBM_workflow.add_nodes([node_3T_MNI_reg_aladin2])

# Node 22
OHBM_workflow.add_nodes([node_3T_MNI_reg_f3d])

# Node 23-1
OHBM_workflow.add_nodes([node_subjMNI_resample])

# Node 23-2
OHBM_workflow.add_nodes([node_1mm_subjMNI_resample])

# Node 24
OHBM_workflow.add_nodes([node_mRs_highres_bet])

# Node 25
OHBM_workflow.add_nodes([node_mRs_despot1])


OHBM_workflow.connect(node_mRs_highres_bet, 'resampled_file', node_reg_aladin, 'floating_img')
OHBM_workflow.connect(node_mRs_highres_bet, 'resampled_file', node_reg_aladin2, 'floating_img')
OHBM_workflow.connect(node_reg_aladin, 'output_affine_transformation', node_reg_aladin2, 'input_affine_trans')

OHBM_workflow.connect(node_mRs_despot1, 'resampled_file', node_reg_aladin, 'reference_img')
OHBM_workflow.connect(node_mRs_despot1, 'resampled_file', node_reg_aladin2, 'reference_img')
OHBM_workflow.connect(node_mRs_despot1, 'resampled_file', node_reg_f3d, 'reference_img')

OHBM_workflow.connect(node_mRs_highres_bet, 'resampled_file', node_reg_f3d, 'floating_img')
OHBM_workflow.connect(node_reg_aladin2, 'output_affine_transformation', node_reg_f3d, 'option_affine')

OHBM_workflow.connect(node_mRs_despot1, 'resampled_file', node_new_reg_transform, 'reference_img')
OHBM_workflow.connect(node_reg_aladin2, 'output_affine_transformation', node_new_reg_transform, 'disp')

OHBM_workflow.connect(node_mRs_despot1, 'resampled_file', node_cpp_reg_transform, 'reference_img')
OHBM_workflow.connect(node_reg_f3d, 'outfile_cpp', node_cpp_reg_transform, 'cpp2def')

OHBM_workflow.connect(node_mRs_despot1, 'resampled_file', node_def_reg_transform, 'reference_img')
OHBM_workflow.connect(node_cpp_reg_transform, 'output_nlin_deformation', node_def_reg_transform, 'def2disp')

OHBM_workflow.connect(node_new_reg_transform, 'output_disp2', node_aff_4d_convert, 'in_nii')
OHBM_workflow.connect(node_def_reg_transform, 'output_nlin_displacement', node_nlin_4d_convert, 'in_nii')

OHBM_workflow.connect(node_nlin_4d_convert, 'output_4d_nii', node_fsl_subtract, 'in_file')
OHBM_workflow.connect(node_aff_4d_convert, 'output_4d_nii', node_fsl_subtract, 'operand_file')

OHBM_workflow.connect(node_fsl_subtract, 'out_file', node_fsl_split, 'in_file')

OHBM_workflow.connect(node_fsl_split, 'out_files', node_select1, 'inlist')
OHBM_workflow.connect(node_select1, 'out', node_fsl_sqr1, 'in_file')

OHBM_workflow.connect(node_fsl_split, 'out_files', node_select2, 'inlist')
OHBM_workflow.connect(node_select2, 'out', node_fsl_sqr2, 'in_file')

OHBM_workflow.connect(node_fsl_split, 'out_files', node_select3, 'inlist')
OHBM_workflow.connect(node_select3, 'out', node_fsl_sqr3, 'in_file')

OHBM_workflow.connect(node_fsl_sqr2, 'out_file', node_merge, 'in1')
OHBM_workflow.connect(node_fsl_sqr3, 'out_file', node_merge, 'in2')

OHBM_workflow.connect(node_fsl_sqr1, 'out_file', node_fsl_add_all, 'in_file') 
OHBM_workflow.connect(node_merge, 'out', node_fsl_add_all, 'operand_files')

OHBM_workflow.connect(node_fsl_add_all, 'out_file', node_fsl_sqrt, 'in_file')

OHBM_workflow.connect(node_mRs_despot1, 'resampled_file', node_3T_MNI_reg_aladin1, 'floating_img')
OHBM_workflow.connect(node_mRs_despot1, 'resampled_file', node_3T_MNI_reg_aladin2, 'floating_img')
OHBM_workflow.connect(node_mRs_despot1, 'resampled_file', node_3T_MNI_reg_f3d, 'floating_img')
OHBM_workflow.connect(node_3T_MNI_reg_aladin2, 'output_affine_transformation', node_3T_MNI_reg_f3d, 'option_affine')

OHBM_workflow.connect(node_3T_MNI_reg_aladin1, 'output_affine_transformation', node_3T_MNI_reg_aladin2, 'input_affine_trans')

OHBM_workflow.connect(node_fsl_sqrt, 'out_file', node_subjMNI_resample, 'floating_img')
OHBM_workflow.connect(node_3T_MNI_reg_f3d, 'outfile_cpp', node_subjMNI_resample, 'cpp_img')

OHBM_workflow.connect(node_fsl_sqrt, 'out_file', node_1mm_subjMNI_resample, 'floating_img')
OHBM_workflow.connect(node_3T_MNI_reg_f3d, 'outfile_cpp', node_1mm_subjMNI_resample, 'cpp_img')

##############
# Draw graph #
##############

#simple
OHBM_workflow.write_graph()

OHBM_workflow.run()
