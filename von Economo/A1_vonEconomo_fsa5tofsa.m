% Main script for generating Economo-Koskinas cortical parcellation
% annotation files and resampling them from fsaverage5 to fsaverage
% using FreeSurfer
%
% Reference:
% Lee H, ..., Huang S, Communications Biology, 2026
%
% Requirements:
% - FreeSurfer installed and configured
% - FreeSurfer fsaverage and fsaverage5 surfaces available
% - ENIGMA Toolbox parcellation CSV file (https://enigma-toolbox.readthedocs.io/en/latest/pages/01.install/index.html)
% - FreeSurfer MATLAB write_annotation function available
%
% Author: Hansol Lee (0000-0003-2112-1197)
% Athinoula A. Martinos Center for Biomedical Imaging,
% Massachusetts General Hospital / Harvard Medical School

clear all; close all; clc;

%% ------------------------------------------------------------
% Settings
% ------------------------------------------------------------

% 
csv_path = '/autofs/cluster/connectome2/Bay8_C2/bids/derivatives/SANDI/ENIGMA/matlab/shared/parcellations/economo_koskinas_fsa5.csv';

%% ------------------------------------------------------------
% Load CSV labels
% ------------------------------------------------------------
labels_all = csvread(csv_path);
nverts_hemi = 10242;

labels_lh = labels_all(1:nverts_hemi);
labels_rh = labels_all(nverts_hemi+1:end);

% Vertex indices (FreeSurfer 0-based index)
vertices_lh = (0:nverts_hemi-1)';
vertices_rh = (0:nverts_hemi-1)';

%% ------------------------------------------------------------
% Build colortable (6 categories)
% ------------------------------------------------------------
ct.numEntries = 6;
ct.orig_tab = 'custom';
ct.struct_names = {'unknown';'Agranular';'FrontalGranular';'Parietal';'Polar';'Granular'};
ct.table = [ ...
    25   25   25   0   25 + 25*2^8 + 25*2^16;      % 0 unknown
    255   0    0   0   255 + 0*2^8 + 0*2^16;       % 1 Agranular
    255 150    0   0   255 + 150*2^8 + 0*2^16;     % 2 Frontal granular
    0   255    0   0   0 + 255*2^8 + 0*2^16;       % 3 Granular
    0     0  255   0   0 + 0*2^8 + 255*2^16;       % 4 Polar/Limbic
    200 100  255   0   200 + 100*2^8 + 255*2^16];  % 5 Dysgranular

%% Convert label indices to color integer values
labels_lh_rgb = ct.table(labels_lh + 1, 5);
labels_rh_rgb = ct.table(labels_rh + 1, 5);

%% Write annotation files
write_annotation('lh.economo_koskinas_fsa5.annot', vertices_lh, labels_lh_rgb, ct);
write_annotation('rh.economo_koskinas_fsa5.annot', vertices_rh, labels_rh_rgb, ct);

%%

dpFs = fullfile('/autofs/cluster/connectome2/Bay8_C2/bids/derivatives/SANDI/freesurfer');
setenv('SUBJECTS_DIR', dpFs);

system(['mri_surf2surf --srcsubject fsaverage5 ' ...
    '--trgsubject fsaverage ' ...
    '--hemi lh ' ...
    '--sval-annot lh.economo_koskinas_fsa5.annot ' ...
    '--tval lh.economo_koskinas_fsaverabe.annot'])

system(['mri_surf2surf --srcsubject fsaverage5 ' ...
    '--trgsubject fsaverage ' ...
    '--hemi rh ' ...
    '--sval-annot rh.economo_koskinas_fsa5.annot ' ...
    '--tval rh.economo_koskinas_fsaverabe.annot'])