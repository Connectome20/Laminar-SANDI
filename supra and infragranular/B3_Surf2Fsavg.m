% Main script for resampling subject-specific supragranular and
% infragranular SANDI surface overlays to the fsaverage surface
% using FreeSurfer mri_surf2surf
%
% Reference:
% Lee H, ..., Huang S, Communications Biology, 2026
%
% Requirements:
% - FreeSurfer installed and configured
% - FreeSurfer recon-all completed for each subject
% - Subject-specific laminar SANDI surface overlays for supragranular and infragranular layers
%
% Author: Hansol Lee (0000-0003-2112-1197)
% Athinoula A. Martinos Center for Biomedical Imaging,
% Massachusetts General Hospital / Harvard Medical School

clear all; close all; clc;

%%

fssubj = 'fs';

sub = {'sub-002'};

for sn = 1

    dproot = '/autofs/cluster/connectome2/Bay8_C2/bids';
    dpFs = fullfile(dproot, sub{sn}, 'anat_process');
    setenv('SUBJECTS_DIR', dpFs);

    % super resolution
    dpSurf = fullfile(dproot, 'derivatives/SANDI/RF_Dis2/granular/supres/pv1/sep', sub{sn});
    dpSurf2Fsavg = fullfile(dproot, 'derivatives/SANDI/RF_Dis2/granular/supres/pv1/surf2Fsavg', sub{sn});
    mkdir(dpSurf2Fsavg)
    
    SANDI_lh_surfs = struct2cell(dir(fullfile(dpSurf,'lh.*.mgz')));
    % SANDI_rh_surfs = struct2cell(dir(fullfile(dpSurf,'rh.*.mgz')));
    
    for ii=1:size(SANDI_lh_surfs,2)
    
        fpSANDI_lh = fullfile(SANDI_lh_surfs{2,ii},SANDI_lh_surfs{1,ii});
        % fpSANDI_rh = fullfile(SANDI_rh_surfs{2,ii},SANDI_rh_surfs{1,ii});

        fpSANDI_lh_n = SANDI_lh_surfs{1,ii};
        ix = strfind(fpSANDI_lh_n, '.');
        sm_n = fpSANDI_lh_n(1:ix(4));

        fpSANDIavg_lh = fullfile(dpSurf2Fsavg, [sm_n 'overlay.mgz']);

        % fpSANDI_rh_n = SANDI_rh_surfs{1,ii};
        % ix = strfind(fpSANDI_rh_n, '.');
        % sm_n = fpSANDI_rh_n(1:ix(4));
        
        % fpSANDIavg_rh = fullfile(dpSurf2Fsavg, [sm_n 'overlay.mgz']);

        cmd = ['mri_surf2surf --srcsubject ' fssubj ' --hemi lh --trgsubject fsaverage --sval ' fpSANDI_lh ' --tval ' fpSANDIavg_lh ' --sfmt paint --cortex'];
        [status, result] = system(cmd, '-echo'); 

        % cmd = ['mri_surf2surf --srcsubject ' fssubj ' --hemi rh --trgsubject fsaverage --sval ' fpSANDI_rh ' --tval ' fpSANDIavg_rh ' --sfmt paint --cortex'];
        % [status, result] = system(cmd, '-echo'); 
    
    end
end
