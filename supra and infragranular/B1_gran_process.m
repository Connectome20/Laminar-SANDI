% Main script for sampling volumetric SANDI maps onto inferred equivolume
% cortical surfaces to obtain supragranular and infragranular layer-specific intensities
%
% Reference:
% Lee H, ..., Huang S, Communications Biology, 2026
%
% Requirements:
% - FreeSurfer installed and configured
% - FreeSurfer recon-all completed for each subject
% - FreeSurfer bbregister registration between T1w and dMRI
% - Inferred equivolume cortical surfaces in subject space (generated from A1_fit_infra_supra_surface.py)
% - Subject-specific SANDI maps in NIfTI format
%
% Author: Hansol Lee (0000-0003-2112-1197)
% Athinoula A. Martinos Center for Biomedical Imaging,
% Massachusetts General Hospital / Harvard Medical School

clear all; close all; clc;

%%

dproot = '/autofs/cluster/connectome2/Bay8_C2/bids';

sub = {'sub-002'};

for sn = 1

    dproot = '/autofs/cluster/connectome2/Bay8_C2/bids';
    dpFs = fullfile(dproot, sub{sn}, 'anat_process');
    setenv('SUBJECTS_DIR', dpFs);
    
    dpdat = fullfile(dproot, 'derivatives/SANDI/RF_Dis2/supres', sub{sn});
    SANDI_maps = struct2cell(dir(fullfile(dpdat,'*SANDI-fit*.nii*')));

    dpgrandat = fullfile(dproot, 'derivatives/SANDI/granular', sub{sn});
    dpout = fullfile(dproot, 'derivatives/SANDI/RF_Dis2/granular', sub{sn});
    mkdir(dpout)

    dpBbr = fullfile(dproot, 'derivatives/bbr', sub{sn});
    fpBbrDat = fullfile(dpBbr, 'diff_bbr_t1w.dat');

    aseg = fullfile(dpFs, 'fs/mri/aseg.mgz');

    isovol_lh_inf = fullfile(dpgrandat, 'isovolume_smooth20/lh.inf.isovolume.global.nocras');
    isovol_rh_inf = fullfile(dpgrandat, 'isovolume_smooth20/rh.inf.isovolume.global.nocras');

    fs_lh_inf001 = fullfile(dpFs, 'fs/surf/lh.inf001');
    fs_lh_inf = fullfile(dpFs, 'fs/surf/lh.inf');
    fs_rh_inf001 = fullfile(dpFs, 'fs/surf/rh.inf001');
    fs_rh_inf = fullfile(dpFs, 'fs/surf/rh.inf');

    % processing inf
    cmd = ['cp ' isovol_lh_inf ' ' fs_lh_inf001];
    [status, result] = system(cmd, '-echo');

    cmd = ['mris_convert --to-scanner ' fs_lh_inf001 ' ' fs_lh_inf001];
    [status, result] = system(cmd, '-echo');

    cmd = ['mris_convert --vol-geom ' aseg ' ' fs_lh_inf001 ' ' fs_lh_inf001];
    [status, result] = system(cmd, '-echo');

    cmd = ['mris_convert --to-tkr ' fs_lh_inf001 ' ' fs_lh_inf001];
    [status, result] = system(cmd, '-echo');

    % cmd = ['cp ' isovol_rh_inf ' ' fs_rh_inf001];
    % [status, result] = system(cmd, '-echo');
    % 
    % cmd = ['mris_convert --to-scanner ' fs_rh_inf001 ' ' fs_rh_inf001];
    % [status, result] = system(cmd, '-echo');
    % 
    % cmd = ['mris_convert --vol-geom ' aseg ' ' fs_rh_inf001 ' ' fs_rh_inf001];
    % [status, result] = system(cmd, '-echo');
    % 
    % cmd = ['mris_convert --to-tkr ' fs_rh_inf001 ' ' fs_rh_inf001];
    % [status, result] = system(cmd, '-echo');

    for ii=1:6

        fpSANDI = fullfile(SANDI_maps{2,ii},SANDI_maps{1,ii});
        vol_n = SANDI_maps{1,ii};
        ix_1 = strfind(vol_n, '_');
        ix_2 = strfind(vol_n, '.');

        if ii == 1
            
            % partial volume fraction
            SANDI_lh_whole_cortex_fraction_maps = fullfile(dpout, 'lh.fractions.whole_cortex.mgz');
            SANDI_lh_inf_sup_fraction_maps = fullfile(dpout, 'lh.fractions.inf_sup.mgz');
        
            SANDI_rh_whole_cortex_fraction_maps = fullfile(dpout, 'rh.fractions.whole_cortex.mgz');
            SANDI_rh_inf_sup_fraction_maps = fullfile(dpout, 'rh.fractions.inf_sup.mgz');
        
            cmd = ['mri_compute_layer_fractions -sdir ' dpFs ' -lh -s fs -nlayers 1 -r 1.5 -FS_names ' fpBbrDat ' ' fpSANDI ' ' SANDI_lh_whole_cortex_fraction_maps];
            [status, result] = system(cmd, '-echo');

            cmd = ['mri_compute_layer_fractions -sdir ' dpFs ' -lh -n inf -s fs -nlayers 2 -r 1.5 -FS_names ' fpBbrDat ' ' fpSANDI ' ' SANDI_lh_inf_sup_fraction_maps];
            [status, result] = system(cmd, '-echo');
        
            % cmd = ['mri_compute_layer_fractions -sdir ' dpFs ' -rh -s fs -nlayers 1 -r 1.5 -FS_names ' fpBbrDat ' ' fpSANDI ' ' SANDI_rh_whole_cortex_fraction_maps];
            % [status, result] = system(cmd, '-echo');
            % 
            % cmd = ['mri_compute_layer_fractions -sdir ' dpFs ' -rh -n inf -s fs -nlayers 2 -r 1.5 -FS_names ' fpBbrDat ' ' fpSANDI ' ' SANDI_rh_inf_sup_fraction_maps];
            % [status, result] = system(cmd, '-echo');
    
        end

        SANDI_lh_whole_cortex_fraction_intensity_maps = fullfile(dpout, ['lh.' vol_n(ix_1(1)+1:ix_2(1)-1) '.frac_intensitiy.whole_cortex.mgz']);
        SANDI_lh_inf_sup_fraction_intensity_maps = fullfile(dpout, ['lh.' vol_n(ix_1(1)+1:ix_2(1)-1) '.frac_intensitiy.inf_sup.mgz']);

        SANDI_rh_whole_cortex_fraction_intensity_maps = fullfile(dpout, ['rh.' vol_n(ix_1(1)+1:ix_2(1)-1) '.frac_intensitiy.whole_cortex.mgz']);
        SANDI_rh_inf_sup_fraction_intensity_maps = fullfile(dpout, ['rh.' vol_n(ix_1(1)+1:ix_2(1)-1) '.frac_intensitiy.inf_sup.mgz']);

        % partial intensity
        cmd = ['mris_compute_layer_intensities -lh -sdir ' dpFs ' -s fs -FS_names -nlayers 1 ' fpSANDI ' ' SANDI_lh_whole_cortex_fraction_maps ' ' fs_lh_inf ' ' SANDI_lh_whole_cortex_fraction_intensity_maps]; 
        [status, result] = system(cmd, '-echo');

        cmd = ['mris_compute_layer_intensities -lh -sdir ' dpFs ' -s fs -FS_names -nlayers 2 ' fpSANDI ' ' SANDI_lh_inf_sup_fraction_maps ' ' fs_lh_inf ' ' SANDI_lh_inf_sup_fraction_intensity_maps];
        [status, result] = system(cmd, '-echo');

        % cmd = ['mris_compute_layer_intensities -rh -sdir ' dpFs ' -s fs -FS_names -nlayers 1 ' fpSANDI ' ' SANDI_rh_whole_cortex_fraction_maps ' ' fs_rh_inf ' ' SANDI_rh_whole_cortex_fraction_intensity_maps]; 
        % [status, result] = system(cmd, '-echo');
        % 
        % cmd = ['mris_compute_layer_intensities -rh -sdir ' dpFs ' -s fs -FS_names -nlayers 2 ' fpSANDI ' ' SANDI_rh_inf_sup_fraction_maps ' ' fs_rh_inf ' ' SANDI_rh_inf_sup_fraction_intensity_maps];
        % [status, result] = system(cmd, '-echo');

    end

    delete(fs_lh_inf001)
    % delete(fs_rh_inf001)

end
