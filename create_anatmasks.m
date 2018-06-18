% This script creates the section of the subject's individual wm anat file
% (i1) and the SPM wm MNI patch (i2): '(i1>0) & (i2>0)'. The same is performed
% for the ventricles. It also combines the subject's individual cortex, wm and
% ventricle segmented anat files into a wholebrain anat file.
% Requires SPM12: www.fil.ion.ucl.ac.uk/spm/software/spm12/
% Path variables below should be set properly.

clear all

%% Basic settings (including path variables)

basename='sald';
rootpath='c:\projects\transfer_learning\SALD_spm\';
path_patchtemplate='c:\projects\matlabtools\custom_spm\MNIpatch\';
fname_liquor_patch='MNI_liquorpatch.nii';
fname_white_patch='MNI_whitepatch.nii';
fls=dir(sprintf('%s0*',rootpath));

%% Initializing SPM

spmfolder = 'c:\projects\matlabtools\spm12';                % Path to SPM.
addpath(spmfolder)
spm('Defaults','FMRI')
spm_jobman('initcfg')

%% Creating WM, CSF, and whole-brain masks

for subnum=1:length(fls)
    subname=fls(subnum).name;
    corename=sprintf('%s_sub%s',basename,subname);
    
    fname_liquor_patch_sub=sprintf('MNI_liquorpatch_%s.nii',corename);
    fname_white_patch_sub=sprintf('MNI_whitepatch_%s.nii',corename);
    anatfolder=sprintf('%s\\%s\\anat\\',rootpath,subname);
    
    copyfile(fullfile(path_patchtemplate,fname_liquor_patch),fullfile(anatfolder,fname_liquor_patch_sub))
    copyfile(fullfile(path_patchtemplate,fname_white_patch),fullfile(anatfolder,fname_white_patch_sub))
    
    cortexfile_name=fullfile(anatfolder,sprintf('wc1%s_anat.nii',corename));
    whitefile_name=fullfile(anatfolder,sprintf('wc2%s_anat.nii',corename));
    ventrfile_name=fullfile(anatfolder,sprintf('wc3%s_anat.nii',corename));
    wholebrain_name=fullfile(anatfolder,sprintf('wwholebrainanat_%s.nii',corename));
    
    
    spm_imcalc({whitefile_name, fullfile(anatfolder,fname_white_patch_sub)},...
        fullfile(anatfolder,fname_white_patch_sub),'(i1>0) & (i2>0)');
    spm_imcalc({ventrfile_name, fullfile(anatfolder,fname_liquor_patch_sub)},...
        fullfile(anatfolder,fname_liquor_patch_sub),'i1>0 & i2>0');   
    spm_imcalc({cortexfile_name,whitefile_name,ventrfile_name},wholebrain_name,'i1|i2|i3');
end


