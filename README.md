# Resolution recovery in routine clinical neuroimaging data

This code enables resolution recovery in clinial grade neuroimaging data. It can process both MRI (denoising or super-resolution) and CT (denoising) data.

**Super-resolution**: The code can reconstruct high-resolution, isotropic data from clinical MR scans (i.e. thick-sliced), of arbitrary orientation and MR contrast. For example, given a T1, a T2 and a FLAIR scan with large slice-thickness, three 1 mm isotropic scans can be recovered. 

**Denoising**: The code can remove noise from MR and CT scans.

The code is based on the method described in this paper:

     Brudfors M, Balbastre Y, Nachev P, Ashburner J.
     MRI Super-Resolution Using Multi-channel Total Variation.
     In Annual Conference on Medical Image Understanding and Analysis
     2018 Jul 9 (pp. 217-228). Springer, Cham.
     
The most up-to-date PDF version of the paper is available from https://arxiv.org/abs/1810.03422.

## Dependencies

This project has strong dependencies on SPM12 and its `Shoot` toolbox. Both of them should be added to Matlab's path. The most recent version of SPM can be downloaded from [www.fil.ion.ucl.ac.uk/spm](http://www.fil.ion.ucl.ac.uk/spm/). If you get error messages when running the code, it is probably because your SPM version is too old.

## Example 1: Super-resolve a set MRIs of one subject

~~~~
% Read some MRI NIfTIs
dir_data = '/pth/to/nii_data';
Nii      = nifti(spm_select('FPList',dir_data,'^.*\.nii$'));

% Super-resolve the MRIs
spm_mtv_preproc('InputImages',Nii,'Method','superres');
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

## Example 3: Denoise a CT image

~~~~
% Read a CT NIfTI
dir_data = '/pth/to/nii_data';
Nii      = nifti(spm_select('FPList',dir_data,'^.*\.nii$'));

% Denoise the CT
spm_mtv_preproc('InputImages',Nii,'Modality','CT');
~~~~

The denoised images will written to the 'out' folder, prefixed 'den_'.
