# MTVprocess3D

Multi-channel total variation (MTV) denoising or super-resolution, of magnetic resonance (MR) images. 

Requires that the SPM software, and our Auxiliary toolbox, is on the MATLAB path. SPM is available from
https://www.fil.ion.ucl.ac.uk/spm/software/spm12/ and the Auxiliary toolbox from https://github.com/WTCN-computational-anatomy-group/auxiliary-functions.

The general principles of the methods are described in the following paper:

     Brudfors M, Balbastre Y, Nachev P, Ashburner J.
     MRI Super-Resolution Using Multi-channel Total Variation.
     In Annual Conference on Medical Image Understanding and Analysis
     2018 Jul 9 (pp. 217-228). Springer, Cham.

## Example 1: Denoising MRIs

~~~~
dir_MTVprocess3D      = '/path/to/MTVprocess3D';
dir_auxiliary_toolbox = '/path/to/auxiliary-functions';

addpath(genpath(dir_MTVprocess3D));
addpath(dir_auxiliary_toolbox);

% Read some MRI NIfTIs
dir_data = '/pth/to/nii_data';
Nii = nifti(spm_select('FPList',dir_data,'^.*\.nii$'));

% Denoise the MRIs
MTVprocess3D('InputImages',Nii);
~~~~
