# Super-resolution and denoising of routine clinical magnetic resonance scans

This code enables reconstructing high-resolution, isotropic data from clinical magnetic resonance (MR) scans, of arbitrary orientation and MR contrast. For example, given a T1, a T2 and a FLAIR scan with large slice-thickness, three 1 mm isotropic scans can be recovered. Furthermore, the method can also be used to denoise noisy MR scans.

The method is described in detail in the following paper:

     Brudfors M, Balbastre Y, Nachev P, Ashburner J.
     MRI Super-Resolution Using Multi-channel Total Variation.
     In Annual Conference on Medical Image Understanding and Analysis
     2018 Jul 9 (pp. 217-228). Springer, Cham.
     
The most up-to-date PDF version of this paper is available from https://arxiv.org/abs/1810.03422.

## Requirements

The code requires that the SPM software is on the MATLAB path. SPM is available from https://www.fil.ion.ucl.ac.uk/spm/software/spm12/.

## Example 1: Super-resolve a set MRIs of one subject

~~~~
% Read some MRI NIfTIs
dir_data = '/pth/to/nii_data';
Nii      = nifti(spm_select('FPList',dir_data,'^.*\.nii$'));

% Super-resolve the MRIs
spm_mtv_preproc('InputImages',Nii,'Method','superres';
~~~~

The super-resolved images will written to the 'out' folder, prefixed 'sr_'.

## Example 2: Denoise a set MRIs of one subject

~~~~
% Read some MRI NIfTIs
dir_data = '/pth/to/nii_data';
Nii      = nifti(spm_select('FPList',dir_data,'^.*\.nii$'));

% Denoise the MRIs
spm_mtv_preproc('InputImages',Nii);
~~~~

The denoised images will written to the 'out' folder, prefixed 'den_'.
