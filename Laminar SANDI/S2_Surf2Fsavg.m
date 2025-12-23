% Main script for resamping subject-specific laminar SANDI surface overlays to the fsaverage surface
% using FreeSurfer mri_surf2surf
%
% Reference:
% Lee H, ...,  Huang S, Communications Biology, 2026
%
% Requirements:
% - FreeSurfer installed and configured
% - FreeSurfer recon-all completed for each subject
% - Subject-specific SANDI surface overlays in FreeSurfer .mgz format
%
% Author: Hansol Lee (0000-0003-2112-1197)
% Athinoula A. Martinos Center for Biomedical Imaging,
% Massachusetts General Hospital / Harvard Medical School

clear all; close all; clc;

%%

dproot = '/autofs/cluster/connectome2/Bay8_C2/bids';
dpsub = fullfile(dproot,'derivatives/SANDI/RF_Dis2/surf/surf2Fsavg');
fssubj = 'fs8';

% Get list of sub folders
dirInfo = dir(dpsub);
subFolders = dirInfo([dirInfo.isdir]);  % keep only directories
subFolders = subFolders(~ismember({subFolders.name}, {'.', '..'}));  % remove . and ..
sub = {subFolders.name};  % extract folder names as cell array

isForce = 1;

for sn = 1:numel(sub)

    subj = sub{sn}

    dpFs = fullfile(dproot, sub{sn}, 'anat_process');
    setenv('SUBJECTS_DIR', dpFs);

    dpSurf = fullfile(dproot, 'derivatives/SANDI/RF_Dis2/supres/surf', sub{sn});
    dpSurf2Fsavg = fullfile(dproot, 'derivatives/SANDI/RF_Dis2/supres/surf/surf2Fsavg', sub{sn});

    lh_list = dir(fullfile(dpSurf,'lh.*.mgz'));
    rh_list = dir(fullfile(dpSurf,'rh.*.mgz'));
    
    for ii=1:numel(lh_list)
    
        fp_lh = fullfile(lh_list(ii).folder, lh_list(ii).name);
        name_lh = lh_list(ii).name;
        ix_lh = strfind(name_lh, '.');
        mapname = erase(name_lh(1:ix_lh(3)-1), 'lh.');
        out_lh = fullfile(dpSurf2Fsavg, ['lh.' mapname '.fsavg.overlay.mgz']);

        name_rh = rh_list(ii).name;
        fp_rh = fullfile(rh_list(ii).folder, name_rh);
        out_rh = fullfile(dpSurf2Fsavg, ['rh.' mapname '.fsavg.overlay.mgz']);

        if ~exist(out_lh, 'file') || isForce
            if ~exist(dpSurf2Fsavg, 'dir')
                mkdir(dpSurf2Fsavg)
            end

            cmd = ['mri_surf2surf --srcsubject ' fssubj ' --hemi lh --trgsubject fsaverage --sval ' fp_lh ' --tval ' out_lh ' --sfmt paint --cortex'];
            [status, result] = system(cmd, '-echo'); 

            cmd = ['mri_surf2surf --srcsubject ' fssubj ' --hemi rh --trgsubject fsaverage --sval ' fp_rh ' --tval ' out_rh ' --sfmt paint --cortex'];
            [status, result] = system(cmd, '-echo'); 

        else

            fprintf('Skipping fsaverage ', mapname);

        end

    end
end