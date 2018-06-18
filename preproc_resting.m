% This script regresses out head-motion paramters, the mean WM, CSF, and
% whole-brain signals, and band-pass filters the residual time courses from
% all GM voxels using a combination of temporal high-pass (based on the
% regression of ninth-order discrete cosine transform basis set) and low-pass
% (bidirectional 12th-order Butterworth IIR) filters to retain signals only
% within the range of 0.009 and 0.08 Hz.
% Requires SPM12: www.fil.ion.ucl.ac.uk/spm/software/spm12/
% and Tools for NIfTI and ANALYZE image:
% https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
% Path variables below should be set properly.

clear all

%% Initializing SPM

spmfolder = 'c:\projects\matlabtools\spm12';                % Path to SPM.
addpath(spmfolder)
addpath([spm('Dir'),'\external\fieldtrip\preproc\']);
spm('Defaults','FMRI')
spm_jobman('initcfg')

%% Basic settings (including path variables)

basename = 'sald';
rootpath = 'c:\projects\transfer_learning\SALD_spm\';       % Path to preprocessed data.
fls = dir(sprintf('%s0*',rootpath));

%% Preprocessing parameters

fwhm = 8;
if fwhm>0
    topreproc_prefix=sprintf('s%gwcr',fwhm);
else
    topreproc_prefix='r';
end

% Regression params

reg_motion = 1;
reg_white = 1;
reg_ventric = 1;
reg_global = 1;

% Temporal filter params

filter= 1; 
filt_order = 12;
Fbp = [0.009 0.08];

%% Preprocessing

for subnum=1:length(fls)
    
    % Subject-specific directories
    
    subname=fls(subnum).name;
    disp(sprintf('%s',subname));
    processfolder=sprintf('%s\\%s\\process\\',rootpath,subname);
    anatfolder=sprintf('%s\\%s\\anat\\',rootpath,subname);
    restingfolder=sprintf('%s\\%s\\process_resting\\',rootpath,subname);
    if ~isdir(restingfolder)
        mkdir(restingfolder)
    end
    corename=sprintf('%s_sub%s',basename,subname);
    
    % Loading paramfile
    
    paramfile=fullfile(processfolder,'nifti_params.mat'); %'TR','dynums','slinum','ht','rawfnames'
    load(paramfile)
    disp(sprintf('TR set to %g',TR))
    Fs=1/TR;
    act_length=dynums;
    
    % Meanfile and anatomical masks
    
    meanfmri_file=sprintf('%swcmean%s_rest.nii',processfolder,corename);
    funcspace=spm_vol(meanfmri_file);
    M=funcspace.mat;
    
    ventrfile_patch_name=sprintf('%swMNI_liquorpatch_%s.nii',anatfolder,corename);
    whitefile_patch_name=sprintf('%swMNI_whitepatch_%s.nii',anatfolder,corename);
    cortexfile_name=sprintf('%sthr0.1_wc1%s_anat.nii',anatfolder,corename);
    
    wholebrain_name=sprintf('%swwholebrainanat_%s.nii',anatfolder,corename);
    
    whitematter_vol=spm_vol(whitefile_patch_name);
    exp_mask_vol=spm_vol(cortexfile_name);
    ventric_vol=spm_vol(ventrfile_patch_name);
    wholebrain_vol=spm_vol(wholebrain_name);
    
    % Loading motion params to create motion regressors
    
    if reg_motion
        realignparamfile=sprintf('%s\\rp_%s_rest.txt',processfolder,corename);
        regressors=load(realignparamfile);
        regressors=regressors(1:act_length,:);  
    else
        regressors=[];
    end
  
     
    % Loading functional images and creating a mask
    
    funcfile=sprintf('%s%s%s_rest.nii',processfolder,topreproc_prefix,corename);
    funcimg=load_untouch_nii(funcfile);
    funci=funcimg.img;
    
    vXYZi=[1:prod(funcimg.hdr.dime.dim(2:4))]';
    [x, y, z]=ind2sub(funcimg.hdr.dime.dim(2:4),vXYZi);
    vXYZ=[x y z]';
    funci=permute(funci,[4 1 2 3]);
    Y=single(reshape(funci,size(funci,1),[]));
    Y=Y(1:act_length,:);
    
    mask=find(~any(isnan(Y),1) & all(Y,1));
    Y=Y(:,mask);
    vXYZi=vXYZi(mask);
    vXYZ=vXYZ(:,mask);
    
    j = exp_mask_vol.mat\M*[vXYZ;ones(1,length(vXYZ))]; % Coordinates in mask image
    
    % Loading white matter mask to create a white matter signal regressor
    
    if reg_white
        white_matter = spm_get_data(whitematter_vol,j,false) > 0;
        regressors=[regressors mean(Y(:,white_matter),2)];
    end
    
    % Loading ventricle mask to create a ventricle signal regressor
    
    if reg_ventric
        ventric_data = spm_get_data(ventric_vol,j,false) > 0;
        regressors=[regressors mean(Y(:,ventric_data),2)];
    end
    
    % Loading whole brain mask to create global signal regressor
    
    if reg_global
        global_data = spm_get_data(wholebrain_vol,j,false) > 0;
        regressors=[regressors mean(Y(:,global_data),2)];
    end
    
    % Cleaning/extracting cortex mask
    
    cortex_mask = spm_get_data(exp_mask_vol,j,false) > 0;
    Y=double(Y(:,cortex_mask));
    vXYZ=vXYZ(:,cortex_mask);
    vXYZi=vXYZi(cortex_mask);
    
    % Regressing out artefacts
    
    if ~isempty(regressors)
        fprintf('regressing out artefacts \n')
        regressors=[regressors ones(size(regressors,1),1)];
        regressors_struct = spm_sp('Set',regressors);
        Y=spm_sp('res',regressors_struct,Y);
    end
    
    % Temporal filtering
    
    if filter
        disp('filter start')
        hpf=round(1/Fbp(1));
        k= size(Y,1);
        n       = fix(2*(k*TR)/hpf + 1);
        X0      = spm_dctmtx(k,n);
        X0 = X0(:,2:end);
        Y = Y - X0*(X0'*Y);
        Y=ft_preproc_lowpassfilter(Y',Fs,Fbp(2),filt_order,[],[]);
        Y=Y';
        Y=Y(filt_order:end,:);
    else
        disp('filter not used')
    end
    
    % Saving functional images and coordinates in a matrix form after regressing out artefacts
    
    disp('saving preprocessed functionals')
    outfname=sprintf('%s%s%s_mot%g_white%g_ventric%g_global%g_filt%g.nii'...
        ,restingfolder,topreproc_prefix,corename,reg_motion,reg_white,reg_ventric,reg_global,filter);
    save([outfname '_func.mat'],'Y');
    save([outfname '_coor.mat'],'vXYZ');
    
end