####################
# Import libraries #
####################

import os
import sys
import numpy as np
import nibabel as nib

######################################
## Program Instructions + Arguments ##
######################################

try:
    ## Set input and output
    infile_path = sys.argv[1]
    outfile_path = sys.argv[2]

    print ""
    print ""
    print "-----------------------------"
    print "Nifti Squeezer (Singleton Dimension Removal)"
    print "Converts FLOAT32(5d image) to 4D nifti by removing singleton dimensions"
    print "Based on original vasst-dev tool: convertNiftiFLOAT32to4D.m"
    print "contact: parkpatrickj@gmail.com"
    print "-----------------------------"
    print ""

except:
    print ""
    print ""
    print "-----------------------------"
    print "-----------------------------"
    print "Nifti Squeezer (Singleton Dimension Removal)"
    print "Converts FLOAT32(5d image) to 4D nifti by removing singleton dimensions"
    print "Based on original vasst-dev tool: convertNiftiFLOAT32to4D.m"
    print "contact: parkpatrickj@gmail.com"
    print ""
    print "Usage: python nii_squeezer.py [input] [output]"
    print "e.g. python nii_squeezer.py ~/Documents/affine_displacement_field.nii.gz ~/Documents/affine_displacement_field_4d.nii.gz"
    print "-----------------------------"
    print "-----------------------------"
    print ""
    print ""
    print "Error:"
    print ""

## Load input image
input_img = nib.load(infile_path)

## Get image data as numpy array
data_input_img = input_img.get_data()

## Check old dimensions check new dimensions of data
old_shape = data_input_img.shape
print "input matrix shape: ", old_shape

#read in old header
old_header = input_img.header

## Squeeze data
squeezed_input = np.squeeze(data_input_img)

## Check new dimensions of data
new_shape = squeezed_input.shape
print "output matrix shape: ", new_shape

## Create Nifti image from squeezed data and old header
no_aff_img = nib.nifti1.Nifti1Image(squeezed_input, None, header=old_header)


## save new image
nib.save(no_aff_img, outfile_path)



