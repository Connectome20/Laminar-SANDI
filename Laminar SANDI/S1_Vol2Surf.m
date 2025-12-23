% Main script for laminar projection of SANDI maps onto cortical surfaces
% using FreeSurfer mri_vol2surf with 5% cortical depth increments 
%
% Reference:
% Lee H, ..., Huang S, Communications Biology, 2026
%
% Requirements:
% - FreeSurfer installed and configured
% - FreeSurfer recon-all completed
% - Valid FreeSurfer bbregister registration between T1W and dMRI
% - SANDI fitting results in NIfTI format
%
% Author: Hansol Lee (0000-0003-2112-1197)
% Athinoula A. Martinos Center for Biomedical Imaging,
% Massachusetts General Hospital / Harvard Medical School

clear all; close all; clc;

%%

dproot = '/autofs/cluster/connectome2/Bay8_C2/bids';
dpsub = fullfile(dproot,'derivatives/processed_dwi');
fssubj = 'fs8';

dirInfo = dir(dpsub);
subFolders = dirInfo([dirInfo.isdir]);  % keep only directories
subFolders = subFolders(~ismember({subFolders.name}, {'.', '..'}));  % remove . and ..
sub = {subFolders.name};  % extract folder names as cell array

isForce = 0;

for sn = 1:numel(sub)

    subj = sub{sn}

    dpFs = fullfile(dproot, sub{sn}, 'anat_process');
    setenv('SUBJECTS_DIR', dpFs);

    dpdat = fullfile(dproot, 'derivatives/SANDI/RF_Dis2/supres', sub{sn});
    dpBbr = fullfile('/autofs/space/rhapsody_001/users/Bay8_C2/bids/derivatives/supres/bbr', sub{sn});

    dpSurf = fullfile(dproot, 'derivatives/SANDI/RF_Dis2/supres/surf', sub{sn});

    fpBbrDat = fullfile(dpBbr, 'diff_bbr_t1w.dat');

    SANDI_maps = struct2cell(dir(fullfile(dpdat,'*SANDI-fit*.nii*')));

    if ~isfolder(fullfile(dpFs,fssubj)) || ~isfile(fpBbrDat) || isempty(SANDI_maps)
        fprintf('Missing data for %s\n', subj);
        continue
    end

    threshold = 0:0.05:1;

    for ii=1:6

        fpSANDI = fullfile(SANDI_maps{2,ii},SANDI_maps{1,ii});
        vol_n = SANDI_maps{1,ii};
        ix_1 = strfind(vol_n, '_');
        ix_2 = strfind(vol_n, '.');

        mapname = vol_n(ix_1(1)+1:ix_2(1)-1);

        for ndep = 1:21

            percent = int2str(100 - threshold(ndep) * 100);
    
            fpSANDILoverlay = fullfile(dpSurf, ['lh.' mapname '.' percent '.overlay.mgz']);
            fpSANDIRoverlay = fullfile(dpSurf, ['rh.' mapname '.' percent '.overlay.mgz']);
            
            if ~exist(fpSANDILoverlay, 'file') || ~exist(fpSANDIRoverlay, 'file') || isForce
                if ~exist(dpSurf, 'dir')
                    mkdir(dpSurf)
                end
                
                fprintf('Running projection for %s at %.2f\n', vol_n, threshold(ndep));

                cmd = ['mri_vol2surf --srcsubject ' fssubj ' --hemi lh --projfrac ' num2str(threshold(ndep)) ' --surf white --srcreg ' fpBbrDat ' --src ' fpSANDI ' --out ' fpSANDILoverlay ' --cortex'];
                [status, result] = system(cmd, '-echo'); 
        
                cmd = ['mri_vol2surf --srcsubject ' fssubj ' --hemi rh --projfrac ' num2str(threshold(ndep)) ' --surf white --srcreg ' fpBbrDat ' --src ' fpSANDI ' --out ' fpSANDIRoverlay ' --cortex'];
                [status, result] = system(cmd, '-echo'); 
            else
                fprintf('Skipping %s (%.2f): overlays already exist\n', vol_n, threshold(ndep));
            end
        end
    end
end