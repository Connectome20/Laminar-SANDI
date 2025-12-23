% Main script for separating supragranular and infragranular SANDI fraction-weighted surface intensity maps
% derived from partial-volume laminar sampling
%
% Reference:
% Lee H, ..., Huang S, Communications Biology, 2026
%
% Requirements:
% - FreeSurfer installed and configured
% - FreeSurfer recon-all completed for each subject
% - Fraction-weighted laminar SANDI intensity maps (generated from
% B1_gran_process.m)
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
    dpSurf = fullfile(dproot, 'derivatives/SANDI/RF_Dis2/granular/supres/pv1', sub{sn});
    dpSurfsep = fullfile(dproot, 'derivatives/SANDI/RF_Dis2/granular/supres/pv1/sep', sub{sn});
    mkdir(dpSurfsep)
    
    SANDI_lh_inf_sup_surfs = struct2cell(dir(fullfile(dpSurf,'lh.*frac_intensity.inf_sup.mgz')));
    SANDI_lh_whole_surfs = struct2cell(dir(fullfile(dpSurf,'lh.*frac_intensity.whole_cortex.mgz')));
    
    for ii=1:size(SANDI_lh_inf_sup_surfs,2)
    
        fpSANDI_lh_inf_sup = fullfile(SANDI_lh_inf_sup_surfs{2,ii},SANDI_lh_inf_sup_surfs{1,ii});
        fpSANDI_lh_whole = fullfile(SANDI_lh_whole_surfs{2,ii},SANDI_lh_whole_surfs{1,ii});
    
        fpSANDI_lh_inf_sup_n = SANDI_lh_inf_sup_surfs{1,ii};
        ix = strfind(fpSANDI_lh_inf_sup_n, '.');
        sm_n = fpSANDI_lh_inf_sup_n(ix(1):ix(2));

        [lh_inf_sup, lh_M, lh_mr_parms] = load_mgh(fpSANDI_lh_inf_sup);
        [lh_whole, ~, ~] = load_mgh(fpSANDI_lh_whole);

        lh_inf = lh_inf_sup(:,:,:,2);
        lh_sup = lh_inf_sup(:,:,:,3);
        lh_whole = lh_whole(:,:,:,2);

        lh_inf_out = fullfile(dpSurfsep, ['lh' sm_n 'frac_intensity.inf.mgz']);
        save_mgh(lh_inf, lh_inf_out, lh_M, lh_mr_parms);

        lh_sup_out = fullfile(dpSurfsep, ['lh' sm_n 'frac_intensity.sup.mgz']);
        save_mgh(lh_sup, lh_sup_out, lh_M, lh_mr_parms);

        lh_whole_out = fullfile(dpSurfsep, ['lh' sm_n 'frac_intensity.whole.mgz']);
        save_mgh(lh_whole, lh_whole_out, lh_M, lh_mr_parms);

    end
end