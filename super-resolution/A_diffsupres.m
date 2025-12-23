% Main script for super-resolution reconstruction of diffusion MRI data using T1w guidance
%
% Reference:
% Lee H, ..., Huang S, Communications Biology, 2026
%
% Requirements:
% - FreeSurfer installed and configured
% - FreeSurfer recon-all completed for each subject
% - FSL (flirt, fslmaths, convert_xfm)
% - MRtrix3 (mrconvert)
% - mrsupres_voxel (MRtrix3-based super-resolution tool; https://github.com/yixinma9/super_resolution)
% - Preprocessed diffusion MRI data in NIfTI format
% - Brain mask in diffusion space
%
% Author: Hansol Lee (0000-0003-2112-1197)
% Athinoula A. Martinos Center for Biomedical Imaging,
% Massachusetts General Hospital / Harvard Medical School

clear all; close all; clc;

%%

delete(gcp('nocreate'))
my_pool = parpool(24);

rf = '/autofs/cluster/connectome2/Bay8_C2/bids';
dpsub = fullfile(rf,'derivatives/processed_dwi');

% Get list of sub folders
dirInfo = dir(dpsub);
subFolders = dirInfo([dirInfo.isdir]);  % keep only directories
subFolders = subFolders(~ismember({subFolders.name}, {'.', '..'}));  % remove . and ..
sub = {subFolders.name};  % extract folder names as cell array

isForce = 0;

for ii = 1:numel(sub)

dproot = '/autofs/cluster/connectome2/Bay8_C2/bids';
dpFs = fullfile(dproot, sub{ii}, 'anat_process');
setenv('SUBJECTS_DIR', dpFs);

dpout = '/autofs/space/rhapsody_001/users/Bay8_C2/bids';

dpdata = fullfile(dproot,'derivatives/processed_dwi',sub{ii});
dpmask = fullfile(dproot,'derivatives/SANDI/dilated_mask',sub{ii});
dpsupres = fullfile(dpout, 'derivatives/supres', sub{ii});
dpbbr = fullfile(dpout, 'derivatives/supres/bbr', sub{ii});
dpdwi_temp = fullfile(dpout, 'derivatives/supres', sub{ii}, 'temp_dwi');
dpdwisupres_temp = fullfile(dpout, 'derivatives/supres', sub{ii}, 'temp_dwi_supres');

% Define required input file paths
fpT1wMgz = fullfile(dpFs, 'fs/mri', 'brain.mgz');
fpmask = fullfile(dpmask, [sub{ii} '_brain_mask_dil.nii.gz']);
fpS0s  = fullfile(dproot, sub{ii}, 'dwi_process/s13_SANDI/b0s_PA.nii.gz');
fpdiff = fullfile(dpdata, [sub{ii} '_dwi.nii.gz']);

% Check if all required inputs exist
input_exists = exist(fpT1wMgz, 'file') && ...
               exist(fpmask, 'file') && ...
               exist(fpS0s,  'file') && ...
               exist(fpdiff, 'file');

output_exists = ~isempty(dir(fullfile(dpsupres, [sub{ii} '_dwi_upsample.nii.gz'])));

if ~input_exists
    fprintf('Skipping %s (missing input files)\n', sub{ii});
    continue;
end

if output_exists && ~isForce
    fprintf('Skipping %s (already processed)\n', sub{ii});
    continue;
end
stop
mkdir(dpsupres)
mkdir(dpbbr)
mkdir(dpdwi_temp)
mkdir(dpdwisupres_temp)

% mgz to nii
fpT1wMgz = fullfile(dpFs, 'fs/mri', 'brain.mgz');
fpT1wNii = fullfile(dpbbr, 't1w_fs.nii.gz');
cmd = ['mri_convert ' fpT1wMgz ' ' fpT1wNii ' -odt float'];
[status, result] = system(cmd, '-echo');

% upsample mask
fpmask = fullfile(dpmask, [sub{ii} '_brain_mask_dil.nii.gz']);
fpupmask = fullfile(dpsupres, [sub{ii} '_brain_mask.nii.gz']);
cmd = ['mri_convert --upsample 2 ' fpmask ' ' fpupmask ' -rt nearest'];
[status, result] = system(cmd, '-echo');

% averaging S0
fpS0s = fullfile(dproot, sub{ii}, 'dwi_process/s13_SANDI/b0s_PA.nii.gz');
fpS0 = fullfile(dproot, sub{ii}, 'dwi_process/s13_SANDI/b0_avg.nii.gz');
cmd = ['fslmaths ' fpS0s ' -Tmean ' fpS0];
[status, result] = system(cmd, '-echo'); 

% upsample s0
fpS0 = fullfile(dproot, sub{ii}, 'dwi_process/s13_SANDI/b0_avg.nii.gz');
fpupS0 = fullfile(dpsupres, [sub{ii} '_b0_avg.nii.gz']);
cmd = ['mri_convert --upsample 2 ' fpS0 ' ' fpupS0];        
[status, result] = system(cmd, '-echo');

% bbr
fpBbrDat = fullfile(dpbbr, 'diff_bbr_t1w.dat');
fpBbrLta = fullfile(dpbbr, 'diff_bbr_t1w.lta');
fpBbrMat = fullfile(dpbbr, 'diff_bbr_t1w.mat');
cmd = ['bbregister --s fs --mov ' fpupS0 ' --reg ' fpBbrDat ' --lta ' fpBbrLta ' --fslmat ' fpBbrMat ' --dti --6'];
[status, result] = system(cmd, '-echo');

% diff to t1w
fpDiff2T1w = fullfile(dpbbr, 'diff_bbr_t1w.nii.gz');   
cmd = ['flirt -in ' fpupS0 ' -ref ' fpT1wNii ' -applyxfm -init ' fpBbrMat ' -out ' fpDiff2T1w ' -interp spline']; 
[status, result] = system(cmd, '-echo');

% inv aff
fpBbrMatInv = fullfile(dpbbr, 'diff_bbr_t1w_inv.mat');
cmd = ['convert_xfm -omat ' fpBbrMatInv ' -inverse ' fpBbrMat];
[status, result] = system(cmd, '-echo'); 

% t1w to diff
fpT1w2Diff = fullfile(dpbbr, 'diff_bbr_t1w_inv.nii.gz');   
cmd = ['flirt -in ' fpT1wNii ' -ref ' fpupS0 ' -applyxfm -init ' fpBbrMatInv ' -out ' fpT1w2Diff ' -interp spline']; 
[status, result] = system(cmd, '-echo'); 

% convert mask and T1 reference to datatype of float64
fpT1w2Diff_float64 = fullfile(dpbbr, 'diff_bbr_t1w_inv_float64.nii.gz');   
cmd = ['mrconvert ' fpT1w2Diff ' ' fpT1w2Diff_float64 ' -datatype float64 -force']; 
[status, result] = system(cmd, '-echo');

fpupmask_float64 = fullfile(dpsupres, [sub{ii} '_brain_mask_float64.nii.gz']);
cmd = ['mrconvert ' fpupmask ' ' fpupmask_float64 ' -datatype float64 -force']; 
[status, result] = system(cmd, '-echo');

% extract 3D volumes from 4D DWI
fpdiff = fullfile(dpdata, [sub{ii} '_dwi.nii.gz']);
parfor i=0:895
fpdiff_3D = fullfile(dpdwi_temp, [sub{ii} '_dwi_3D_' num2str(i) '.nii.gz']);
cmd = ['mrconvert ' fpdiff ' ' fpdiff_3D ' -coord 3 ' num2str(i) ' -force']; 
[status, result] = system(cmd, '-echo');
end

% mrsupres
parfor i=0:895

    fpdiff_3D = fullfile(dpdwi_temp, [sub{ii} '_dwi_3D_' num2str(i) '.nii.gz']);
    fpdiff_up3D = fullfile(dpdwisupres_temp, [sub{ii} '_dwi_3D_' num2str(i) '.nii.gz']);

    if size(dir(fpdiff_up3D),1) == 0

        mkdir([pwd '/' sub{ii} '/' num2str(i)])

        cd([pwd '/' sub{ii} '/' num2str(i)])

        cmd = ['/autofs/space/linen_001/users/Yixin/mrtrix3/bin/mrsupres_voxel ' fpdiff_3D ' ' fpdiff_up3D ' ' ...
            '-mask ' fpupmask_float64 ' -hires ' fpT1w2Diff_float64 ' -kernel 2 -iter_params 1,2,4,6,8,16 -nthread 64 -force'];
        [status, result] = system(cmd, '-echo');

        cd('../..')

    end
end

% merge 3D volumes to 4D DWI
fpdiff_up3D_all = fullfile(dpdwisupres_temp, [sub{ii} '_dwi_3D_[].nii.gz']);
fpdiff_up4D = fullfile(dpsupres, [sub{ii} '_dwi_upsample.nii.gz']);
cmd = ['mrconvert ' fpdiff_up3D_all ' ' fpdiff_up4D];
[status, result] = system(cmd, '-echo');

% remove the temp
rmdir([pwd '/' sub{ii}],'s')
rmdir(dpdwi_temp,'s')
rmdir(dpdwisupres_temp,'s')

end
