% This script thresholds the subject's individual cortex anat file at 0.1.
% Requires SPM12: www.fil.ion.ucl.ac.uk/spm/software/spm12/
% Path variables below should be set properly.

clear all

%% Basic settings (including path variables)

basename='sald';
rootpath='c:\projects\transfer_learning\SALD_spm\';
fls=dir(sprintf('%s0*',rootpath));

%% Initializing SPM

spmfolder = 'c:\projects\matlabtools\spm12';                % Path to SPM.
addpath(spmfolder)
spm('Defaults','FMRI')
spm_jobman('initcfg')

%% Thresholding canat

for subnum=1:length(fls)
    subname=fls(subnum).name;
    corename=sprintf('%s_sub%s',basename,subname);
    anatfolder=sprintf('%s\\%s\\anat\\',rootpath,subname);
    cortexfile=sprintf('wc1%s_anat.nii',corename);
    spm_imcalc(fullfile(anatfolder,cortexfile),fullfile(anatfolder,['thr0.1_' cortexfile]),'i1>0.1')
end