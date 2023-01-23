# M0_to_PD
This toolbox is a generalization of a mrQ section.
Eliminates the bias field found in M0 image using the linear relationship between M0 and another map (typically T1).
Can be used also for weighted images (PDw as M0 and T1w as T1).

The main function is M0_toPD.m

# Inputs:
* output dir - a directory in which the output PD image will be created.
* M0 file - some PD-weighted image containing a field bias you wish to remove.
* BM file - brain mask. 0 outside the brain, 1 inside the brain.
* Qmap file - some image that possesses a linear relationship with the PD-weighted image.
* seg file - a segmentation image of the same size as the previous images, that contains a segmentation of tissue types - typically GM, WM and CSF.
* Qmap factor - the ratio between the values of the M0 image and the Qmap image. For M0 and T1, it is 1.

# Outputs:
* PD file with a fixed bias field.

# Software requirments:
MATLAB - http://www.mathworks.com/products/matlab/
Vistasoft - https://github.com/vistalab/vistasoft
