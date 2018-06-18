# rsfMRI_fconn

This repository contains scripts that implement the preprocessing of resting-state functional magnetic resonance imaging data and
the calculation of region-of-interest based whole-brain functional connectivity. These scripts were used in the manuscript titled 
'Transfer learning improves resting-state functional connectivity pattern analysis using convolutional neural networks' by Vakli, 
Deák-Meszlényi, Hermann, & Vidnyánszky (submitted for peer review).

Requirements:

* MATLAB versions 7.4 (R2007a) to 9.4 (R2018a; The MathWorks Inc., Natick, MA, USA): https://www.mathworks.com/products/
* SPM12: http://www.fil.ion.ucl.ac.uk/spm/software/spm12/
* Tools for NIfTI and ANALYZE image: https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
* Harvard-Oxford Atlas included in FSL, available from: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Atlases

Raw anatomical and functional NIfTI files should be listed as follows:

* sald_sub031521_rest.nii
* sald_sub031521_anat.nii

where 'sald' is the identifier of the dataset (in this case, the Southwest University Adult lifespan Dataset), coded in the variable 
'basename' throught in the scripts; sub031521 is the subject identifier (in this case, subject #31521 in the SALD dataset) listed in 
the variable 'subnames'; 'rest' and 'anat' are the types of the measurements, coded in the variable 'exptype'. Path variables should
be properly adjusted in the scripts.

To perform data preprocessing and functional connectivity calculation, the scripts should be executed in the following order:

1. 'resting_batch.m' implements motion correction, structural-functional coregistration, spatial normalization, and smoothing.
2. 'create_anatmasks.m' creates individual anatomical masks.
3. 'threshold_canat.m' thresholds each subject's cortical mask.
4. 'preproc_resting.m' regresses out nuisance signals and implements temporal band-pass filtering.
5. 'data_prepare_resting.m' extracts the time courses averaged for each region of interest of the Harvard-Oxford Atlas for each subject separately.
6. 'corr_maps.m' claculates the pairwise linear correlation coefficient between each pair of ROIs for each subject separately.

