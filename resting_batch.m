% This script applies the following preprocessing steps to the T2*-weighted
% resting-state functional images: motion correction (set 'do_realign' to
% 1 and 'do_copy' to 1), structural-functional coregistration (set 'do_coreg'
% to 1), spatial normalization to MNI space (set 'do_normalize' to 1 first; 
% then set 'do_norm_write' to 1 after coregistration has been done) spatial 
% smoothing with a full-width half maximum of 8 mm (set 'do_smooth' to 1; 
% fwhm = 8). Set all the values of the variables below accordingly.
% Requires SPM12: www.fil.ion.ucl.ac.uk/spm/software/spm12/

clear all

%% Basic settings (including path variables)

basename = 'sald';
exptype = 'rest';
rootpath = 'c:\projects\transfer_learning\';        % Study path
rawnifti_folder = [rootpath 'SALD_raw_data_nii\'];  % Path to raw nifti files.
fls = dir(sprintf('%s%s*_%s.nii',rawnifti_folder,basename,exptype));

%% Initializing SPM

spmfolder = 'c:\projects\matlabtools\spm12';                % Path to SPM.
addpath(spmfolder)
codefolder = pwd;                                           % Current folder containing the codes.
spm('Defaults','FMRI')
spm_jobman('initcfg')
rmpath('c:\projects\matlabtools\spm12\toolbox\OldNorm');    % Removing folder from SPM path.
rmpath('c:\projects\matlabtools\spm12\toolbox\OldSeg');     % Removing folder from SPM path.

%% Preprocessing

% Preprocessing parameters

do_realign = 1;                                     % Motion correction.
refsli = 1;                                         % Reference slice for motion correction.
do_copy = 1;                                        % Moving motion corrected functional and anatomical images to the SPM folder for each subject.
do_coreg = 1;                                       % Coregistration.
do_normalize = 1;                                   % Calculating spatial normalization parameters.
do_norm_write = 1;                                  % Spatial normalization.
tonormwrite_prefix = 'cr'; 
do_smooth = 1;                                      % Spatial smoothing.
fwhm = 8;                                           % Full-width half maximum.
tosmooth_prefix = 'wcr';

% Performing the preprocessing

for subnum=1:length(fls)

    fname=fls(subnum).name;
    subname=fname(regexp(fname,'\d'));
    spmfolder=sprintf('%sSALD_spm\\%s\\',rootpath,subname);
    processfolder=sprintf('%sprocess\\',spmfolder);
    anatfolder=sprintf('%sSALD_spm\\%s\\anat\\',rootpath,subname);
    anafile=sprintf('%s_sub%s_anat.nii',basename,subname);
    corename=sprintf('%s_sub%s_%s',basename,subname,exptype);
    
    if ~isdir(processfolder)
        mkdir(processfolder)
    end
    if ~isdir(anatfolder)
        mkdir(anatfolder)
    end
    paramfile=fullfile(processfolder,'nifti_params.mat');
    rawfnames=[];
    dynums=[];
    slinum=[];
    if ~exist(paramfile)
        ht_act=spm_vol(fullfile(rawnifti_folder,fname));
        dynums=length(ht_act);
        slinum=ht_act(1).dim(3);
        rawfname=fname(1:end-4);
        ht=ht_act(1);
        TR=ht(1).private.timing.tspace;
        save(paramfile,'dynums','slinum','ht','rawfname','TR')
        clear ht ht_act
    else
        load(paramfile)
    end
    
    if do_realign
        tempbatch=[];
        for i=1:dynums
            tempbatch{1}.spm.spatial.realign.estwrite.data{1,1}{i,1} = ...
                sprintf('%s%s.nii,%g',rawnifti_folder,rawfname,i);
        end
        tempbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
        tempbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
        tempbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
        tempbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
        tempbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 4;
        tempbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
        tempbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = {''};
        tempbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];%all + mean
        tempbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
        tempbatch{1}.spm.spatial.realign.write.roptions.graphics = 1;
        tempbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
        tempbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
        tempbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
        spm_jobman('run',tempbatch)
    end
    
    if do_copy
        tempbatch=[];
        infile=sprintf('%sr%s.nii',rawnifti_folder,rawfname);
        outfile=sprintf('%sr%s.nii',processfolder,corename);
        movefile(infile,outfile)
        infile=sprintf('%s%s.mat',rawnifti_folder,rawfname);
        outfile=sprintf('%s%s.mat',processfolder,corename);
        movefile(infile,outfile)
        infile=sprintf('%s\\rp_%s.txt',rawnifti_folder,rawfname);
        outfile=sprintf('%s\\rp_%s.txt',processfolder,corename);
        movefile(infile,outfile)
        infile=sprintf('%s\\mean%s.nii',rawnifti_folder,rawfname);
        outfile=sprintf('%s\\mean%s.nii',processfolder,corename);
        movefile(infile,outfile)
        infile=sprintf('%s%s',rawnifti_folder,anafile);
        outfile=sprintf('%s%s',anatfolder,anafile);
        copyfile(infile,outfile)
    end
    
    if do_coreg
        tempbatch=[];
        meanfmri_file=sprintf('%s\\mean%s.nii',processfolder,corename);
        tempbatch{1}.spm.spatial.coreg.estimate.ref{1} = fullfile(anatfolder,anafile);
        tempbatch{1}.spm.spatial.coreg.estimate.source{1} = meanfmri_file;
        
        for i=1:dynums
            tempbatch{1}.spm.spatial.coreg.estimate.other{i,1} = ...
                sprintf('%sr%s.nii,%g',processfolder,corename,i);
        end
        tempbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        tempbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        tempbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        tempbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
        spm_jobman('run',tempbatch)
        
        infile=sprintf('%sr%s.nii',processfolder,corename);
        outfile=sprintf('%scr%s.nii',processfolder,corename);
        movefile(infile,outfile)
        infile=sprintf('%sr%s.mat',processfolder,corename);
        outfile=sprintf('%scr%s.mat',processfolder,corename);
        movefile(infile,outfile)
        infile=sprintf('%smean%s.nii',processfolder,corename);
        outfile=sprintf('%scmean%s.nii',processfolder,corename);
        movefile(infile,outfile)
    end
    
    if do_normalize
        tempbatch=[];
        tpmdir=[spm('Dir') '\tpm'];
        tempbatch{1}.spm.spatial.preproc.channel.vols = {fullfile(anatfolder,anafile)};
        tempbatch{1}.spm.spatial.preproc.channel.write = [0 1];
        ngaus  = [2 2 2 3 4 2];
        native = [1 1 1 0 0 0];
        for c = 1:6 % tissue class c
            tempbatch{1}.spm.spatial.preproc.tissue(c).tpm = {
                fullfile(tpmdir, sprintf('TPM.nii,%d', c))};
            tempbatch{1}.spm.spatial.preproc.tissue(c).ngaus = ngaus(c);
            tempbatch{1}.spm.spatial.preproc.tissue(c).native = [native(c) 0];
            if c < 4
                tempbatch{1}.spm.spatial.preproc.tissue(c).warped = [1 1];
            else
                tempbatch{1}.spm.spatial.preproc.tissue(c).warped = [0 0];
            end
        end
        tempbatch{1}.spm.spatial.preproc.warp.mrf = 1;
        tempbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
        tempbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
        tempbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
        tempbatch{1}.spm.spatial.preproc.warp.fwhm  = 0;
        tempbatch{1}.spm.spatial.preproc.warp.samp = 3;
        tempbatch{1}.spm.spatial.preproc.warp.write = [1 1];
        spm_jobman('run',tempbatch)
        cortexfile=['c1' anafile];
        spm_imcalc(fullfile(anatfolder,cortexfile),fullfile(anatfolder,['thr_' cortexfile]),'i1>0.1')
    end
    
    if do_norm_write
        deffield_file=fullfile(anatfolder,['y_' anafile]);
        
        %%anat
        tempbatch=[];
        tempbatch{1}.spm.spatial.normalise.write.subj.resample{1,1} = fullfile(anatfolder,['m' anafile]);
        tempbatch{1}.spm.spatial.normalise.write.subj.def = {deffield_file};
        tempbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
            78 76 85];
        tempbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
        tempbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
        spm_jobman('run',tempbatch)
        
        %%funci
        tempbatch=[];
        for i=1:dynums
            tempbatch{1}.spm.spatial.normalise.write.subj.resample{i,1} = ...
                sprintf('%s\\%s%s.nii,%g',processfolder,tonormwrite_prefix,corename,i);
        end
        tempbatch{1}.spm.spatial.normalise.write.subj.def = {deffield_file};
        tempbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
            78 76 85];
        tempbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
        tempbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
        spm_jobman('run',tempbatch)
        
        %%meanfunci
        tempbatch=[];
        tempbatch{1}.spm.spatial.normalise.write.subj.resample{1,1} = sprintf('%s\\cmean%s.nii',processfolder,corename);
        tempbatch{1}.spm.spatial.normalise.write.subj.def = {deffield_file};
        tempbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
            78 76 85];
        tempbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
        tempbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
        spm_jobman('run',tempbatch)
    end
    
    if do_smooth
        tempbatch=[];
        for i=1:dynums
            tempbatch{1}.spm.spatial.smooth.data{i,1} = ...
                sprintf('%s%s%s.nii,%g',processfolder,tosmooth_prefix,corename,i);
        end
        tempbatch{1}.spm.spatial.smooth.fwhm = [fwhm fwhm fwhm];
        tempbatch{1}.spm.spatial.smooth.dtype = 0;
        tempbatch{1}.spm.spatial.smooth.prefix =sprintf('s%g', fwhm); 
        spm_jobman('run',tempbatch)
    end
    

    try
        clear tempbatch
    catch
    end
    cd(codefolder)
    
end
