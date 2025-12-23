% Main script for computing group-mean fsaverage surface maps of
% fraction-weighted laminar SANDI metrics for supragranular and
% infragranular layers across subjects
%
% Reference:
% Lee H, ..., Huang S, Communications Biology, 2026
%
% Requirements:
% - FreeSurfer installed and configured
% - Subject-specific SANDI surface overlays resampled to fsaverage space in FreeSurfer .mgz format
%
% Author: Hansol Lee (0000-0003-2112-1197)
% Athinoula A. Martinos Center for Biomedical Imaging,
% Massachusetts General Hospital / Harvard Medical School

clear all; close all; clc;

%%

order = [1 11 2 4 5 6 7 8 9 10 12 13 14 15 16 17 18 19 20 21 3];

dproot = '/autofs/cluster/connectome2/Bay8_C2/bids';

sub = {'sub-002','sub-003','sub-004','sub-005','sub-006','sub-007','sub-008','sub-009','sub-012','sub-014', ...
    'sub-015','sub-016','sub-017','sub-020','sub-035','sub-036','sub-037','sub-038b','sub-039','sub-040','sub-042'};

dpSurf2Fsavg_Overall = fullfile(dproot, 'derivatives/SANDI/RF_Dis2/granular/supres/pv1/surf2Fsavg_Overall');
mkdir(dpSurf2Fsavg_Overall)

sn_all = [1:21];

for tt = 1:21

    sn = sn_all(tt);

    lh_De_all = [];
    lh_Din_all = [];
    lh_fec_all = [];
    lh_fin_all = [];
    lh_fis_all = [];
    lh_Rs_all = [];
    
    rh_De_all = [];
    rh_Din_all = [];
    rh_fec_all = [];
    rh_fin_all = [];
    rh_fis_all = [];
    rh_Rs_all = [];

    dpSurf2Fsavg = fullfile(dproot, 'derivatives/SANDI/23123/granular/supres/pv1/surf2Fsavg', sub{sn});
    
    lh_surfs = struct2cell(dir(fullfile(dpSurf2Fsavg,'lh.*.mgz')));

    for sur = 1:size(lh_surfs,2)

        lh_all(:,sur,tt) = load_mgh(fullfile(lh_surfs{2,sur},lh_surfs{1,sur}));

    end
end

% compute mean
lh_mean = mean(lh_all, 3);

% save out results
[~, M, mr_parms]  = load_mgh(fullfile(lh_surfs{2,1},lh_surfs{1,1}));
fpLhDeinfMean = fullfile(dpSurf2Fsavg_Overall, ['lh.De.frac_intensity.inf.fsavg.mean.mgz']);
save_mgh(lh_mean(:,1), fpLhDeinfMean, M, mr_parms);

[~, M, mr_parms]  = load_mgh(fullfile(lh_surfs{2,1},lh_surfs{1,1}));
fpLhDesupMean = fullfile(dpSurf2Fsavg_Overall, ['lh.De.frac_intensity.sup.fsavg.mean.mgz']);
save_mgh(lh_mean(:,2), fpLhDesupMean, M, mr_parms);

[~, M, mr_parms]  = load_mgh(fullfile(lh_surfs{2,1},lh_surfs{1,1}));
fpLhDewholeMean = fullfile(dpSurf2Fsavg_Overall, ['lh.De.frac_intensity.whole.fsavg.mean.mgz']);
save_mgh(lh_mean(:,3), fpLhDewholeMean, M, mr_parms);


[~, M, mr_parms]  = load_mgh(fullfile(lh_surfs{2,1},lh_surfs{1,1}));
fpLhDininfMean = fullfile(dpSurf2Fsavg_Overall, ['lh.Din.frac_intensity.inf.fsavg.mean.mgz']);
save_mgh(lh_mean(:,4), fpLhDininfMean, M, mr_parms);

[~, M, mr_parms]  = load_mgh(fullfile(lh_surfs{2,1},lh_surfs{1,1}));
fpLhDinsupMean = fullfile(dpSurf2Fsavg_Overall, ['lh.Din.frac_intensity.sup.fsavg.mean.mgz']);
save_mgh(lh_mean(:,5), fpLhDinsupMean, M, mr_parms);

[~, M, mr_parms]  = load_mgh(fullfile(lh_surfs{2,1},lh_surfs{1,1}));
fpLhDinwholeMean = fullfile(dpSurf2Fsavg_Overall, ['lh.Din.frac_intensity.whole.fsavg.mean.mgz']);
save_mgh(lh_mean(:,6), fpLhDinwholeMean, M, mr_parms);


[~, M, mr_parms]  = load_mgh(fullfile(lh_surfs{2,1},lh_surfs{1,1}));
fpLhRsomainfMean = fullfile(dpSurf2Fsavg_Overall, ['lh.Rsoma.frac_intensity.inf.fsavg.mean.mgz']);
save_mgh(lh_mean(:,7), fpLhRsomainfMean, M, mr_parms);

[~, M, mr_parms]  = load_mgh(fullfile(lh_surfs{2,1},lh_surfs{1,1}));
fpLhRsomasupMean = fullfile(dpSurf2Fsavg_Overall, ['lh.Rsoma.frac_intensity.sup.fsavg.mean.mgz']);
save_mgh(lh_mean(:,8), fpLhRsomasupMean, M, mr_parms);

[~, M, mr_parms]  = load_mgh(fullfile(lh_surfs{2,1},lh_surfs{1,1}));
fpLhRsomawholeMean = fullfile(dpSurf2Fsavg_Overall, ['lh.Rsoma.frac_intensity.whole.fsavg.mean.mgz']);
save_mgh(lh_mean(:,9), fpLhRsomawholeMean, M, mr_parms);


[~, M, mr_parms]  = load_mgh(fullfile(lh_surfs{2,1},lh_surfs{1,1}));
fpLhfextrainfMean = fullfile(dpSurf2Fsavg_Overall, ['lh.fextra.frac_intensity.inf.fsavg.mean.mgz']);
save_mgh(lh_mean(:,10), fpLhfextrainfMean, M, mr_parms);

[~, M, mr_parms]  = load_mgh(fullfile(lh_surfs{2,1},lh_surfs{1,1}));
fpLhfextrasupMean = fullfile(dpSurf2Fsavg_Overall, ['lh.fextra.frac_intensity.sup.fsavg.mean.mgz']);
save_mgh(lh_mean(:,11), fpLhfextrasupMean, M, mr_parms);

[~, M, mr_parms]  = load_mgh(fullfile(lh_surfs{2,1},lh_surfs{1,1}));
fpLhfextrawholeMean = fullfile(dpSurf2Fsavg_Overall, ['lh.fextra.frac_intensity.whole.fsavg.mean.mgz']);
save_mgh(lh_mean(:,12), fpLhfextrawholeMean, M, mr_parms);


[~, M, mr_parms]  = load_mgh(fullfile(lh_surfs{2,1},lh_surfs{1,1}));
fpLhfneuriteinfMean = fullfile(dpSurf2Fsavg_Overall, ['lh.fneurite.frac_intensity.inf.fsavg.mean.mgz']);
save_mgh(lh_mean(:,13), fpLhfneuriteinfMean, M, mr_parms);

[~, M, mr_parms]  = load_mgh(fullfile(lh_surfs{2,1},lh_surfs{1,1}));
fpLhfneuritesupMean = fullfile(dpSurf2Fsavg_Overall, ['lh.fneurite.frac_intensity.sup.fsavg.mean.mgz']);
save_mgh(lh_mean(:,14), fpLhfneuritesupMean, M, mr_parms);

[~, M, mr_parms]  = load_mgh(fullfile(lh_surfs{2,1},lh_surfs{1,1}));
fpLhfneuritewholeMean = fullfile(dpSurf2Fsavg_Overall, ['lh.fneurite.frac_intensity.whole.fsavg.mean.mgz']);
save_mgh(lh_mean(:,15), fpLhfneuritewholeMean, M, mr_parms);


[~, M, mr_parms]  = load_mgh(fullfile(lh_surfs{2,1},lh_surfs{1,1}));
fpLhfsomainfMean = fullfile(dpSurf2Fsavg_Overall, ['lh.fsoma.frac_intensity.inf.fsavg.mean.mgz']);
save_mgh(lh_mean(:,16), fpLhfsomainfMean, M, mr_parms);

[~, M, mr_parms]  = load_mgh(fullfile(lh_surfs{2,1},lh_surfs{1,1}));
fpLhfsomasupMean = fullfile(dpSurf2Fsavg_Overall, ['lh.fsoma.frac_intensity.sup.fsavg.mean.mgz']);
save_mgh(lh_mean(:,17), fpLhfsomasupMean, M, mr_parms);

[~, M, mr_parms]  = load_mgh(fullfile(lh_surfs{2,1},lh_surfs{1,1}));
fpLhfsomawholeMean = fullfile(dpSurf2Fsavg_Overall, ['lh.fsoma.frac_intensity.whole.fsavg.mean.mgz']);
save_mgh(lh_mean(:,18), fpLhfsomawholeMean, M, mr_parms);
