% Main script for computing group-mean fsaverage surface maps of
% laminar SANDI metrics across subjects
%
% Reference:
% Lee H, ...,  Huang S, Communications Biology, 2026
%
% Requirements:
% - FreeSurfer installed and configured
% - FreeSurfer recon-all completed for each subject
% - Subject-specific SANDI surface overlays resampled to fsaverage space in FreeSurfer .mgz format
%
% Author: Hansol Lee (0000-0003-2112-1197)
% Athinoula A. Martinos Center for Biomedical Imaging,
% Massachusetts General Hospital / Harvard Medical School

clear all; close all; clc;

%%

order = [1 11 2 4 5 6 7 8 9 10 12 13 14 15 16 17 18 19 20 21 3];

dproot = '/autofs/cluster/connectome2/Bay8_C2/bids';
dpMean = fullfile(dproot,'/derivatives/SANDI/RF_Dis2/supres/surf/surf2Fsavg/GroupMean_Depth');
mkdir(dpMean)

sub = {'sub-002','sub-003','sub-004','sub-005','sub-006','sub-007','sub-008','sub-009','sub-012','sub-014', ...
    'sub-015','sub-016','sub-017','sub-020','sub-035','sub-036','sub-037','sub-038b','sub-039','sub-040','sub-042'};

for dep = 1:length(order)

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

    for sn = 1:21

    dpSurf2Fsavg_Overall = fullfile(dproot, 'derivatives/SANDI/RF_Dis2/supres/surf/surf2Fsavg', sub{sn});

    lh_De_surfs = struct2cell(dir(fullfile(dpSurf2Fsavg_Overall,'lh.*De*.mgz')));
    rh_De_surfs = struct2cell(dir(fullfile(dpSurf2Fsavg_Overall,'rh.*De*.mgz')));
    lh_Din_surfs = struct2cell(dir(fullfile(dpSurf2Fsavg_Overall,'lh.*Din*.mgz')));
    rh_Din_surfs = struct2cell(dir(fullfile(dpSurf2Fsavg_Overall,'rh.*Din*.mgz')));
    lh_fec_surfs = struct2cell(dir(fullfile(dpSurf2Fsavg_Overall,'lh.*fextra*.mgz')));
    rh_fec_surfs = struct2cell(dir(fullfile(dpSurf2Fsavg_Overall,'rh.*fextra*.mgz')));
    lh_fin_surfs = struct2cell(dir(fullfile(dpSurf2Fsavg_Overall,'lh.*fneurite*.mgz')));
    rh_fin_surfs = struct2cell(dir(fullfile(dpSurf2Fsavg_Overall,'rh.*fneurite*.mgz')));
    lh_fis_surfs = struct2cell(dir(fullfile(dpSurf2Fsavg_Overall,'lh.*fsoma*.mgz')));
    rh_fis_surfs = struct2cell(dir(fullfile(dpSurf2Fsavg_Overall,'rh.*fsoma*.mgz')));
    lh_Rs_surfs = struct2cell(dir(fullfile(dpSurf2Fsavg_Overall,'lh.*Rsoma*.mgz')));
    rh_Rs_surfs = struct2cell(dir(fullfile(dpSurf2Fsavg_Overall,'rh.*Rsoma*.mgz')));
    
    lh_De = load_mgh(fullfile(lh_De_surfs{2,order(dep)},lh_De_surfs{1,order(dep)}));
    rh_De = load_mgh(fullfile(rh_De_surfs{2,order(dep)},rh_De_surfs{1,order(dep)}));
    lh_Din = load_mgh(fullfile(lh_Din_surfs{2,order(dep)},lh_Din_surfs{1,order(dep)}));
    rh_Din = load_mgh(fullfile(rh_Din_surfs{2,order(dep)},rh_Din_surfs{1,order(dep)}));
    lh_fec = load_mgh(fullfile(lh_fec_surfs{2,order(dep)},lh_fec_surfs{1,order(dep)}));
    rh_fec = load_mgh(fullfile(rh_fec_surfs{2,order(dep)},rh_fec_surfs{1,order(dep)}));
    lh_fin = load_mgh(fullfile(lh_fin_surfs{2,order(dep)},lh_fin_surfs{1,order(dep)}));
    rh_fin = load_mgh(fullfile(rh_fin_surfs{2,order(dep)},rh_fin_surfs{1,order(dep)}));
    lh_fis = load_mgh(fullfile(lh_fis_surfs{2,order(dep)},lh_fis_surfs{1,order(dep)}));
    rh_fis = load_mgh(fullfile(rh_fis_surfs{2,order(dep)},rh_fis_surfs{1,order(dep)}));
    lh_Rs = load_mgh(fullfile(lh_Rs_surfs{2,order(dep)},lh_Rs_surfs{1,order(dep)}));
    rh_Rs = load_mgh(fullfile(rh_Rs_surfs{2,order(dep)},rh_Rs_surfs{1,order(dep)}));

    lh_De_all = [lh_De_all, lh_De];
    rh_De_all = [rh_De_all, rh_De];
    lh_Din_all = [lh_Din_all, lh_Din];
    rh_Din_all = [rh_Din_all, rh_Din];
    lh_fec_all = [lh_fec_all, lh_fec];
    rh_fec_all = [rh_fec_all, rh_fec];
    lh_fin_all = [lh_fin_all, lh_fin];
    rh_fin_all = [rh_fin_all, rh_fin];
    lh_fis_all = [lh_fis_all, lh_fis];
    rh_fis_all = [rh_fis_all, rh_fis];
    lh_Rs_all = [lh_Rs_all, lh_Rs];
    rh_Rs_all = [rh_Rs_all, rh_Rs];

    end

    % compute mean
    lh_De_mean = mean(lh_De_all, 2);
    rh_De_mean = mean(rh_De_all, 2);
    lh_Din_mean = mean(lh_Din_all, 2);
    rh_Din_mean = mean(rh_Din_all, 2);
    lh_fec_mean = mean(lh_fec_all, 2);
    rh_fec_mean = mean(rh_fec_all, 2);
    lh_fin_mean = mean(lh_fin_all, 2);
    rh_fin_mean = mean(rh_fin_all, 2);
    lh_fis_mean = mean(lh_fis_all, 2);
    rh_fis_mean = mean(rh_fis_all, 2);
    lh_Rs_mean = mean(lh_Rs_all, 2);
    rh_Rs_mean = mean(rh_Rs_all, 2);
    
    % save out results
    [~, M, mr_parms]  = load_mgh(fullfile(lh_De_surfs{2,1},lh_De_surfs{1,1}));
    fpLhDeMean = fullfile(dpMean, ['lh.De.' num2str(str2double(extract(lh_De_surfs{1,order(dep)}, digitsPattern))) '.fsavg.mean.mgz']);
    save_mgh(lh_De_mean, fpLhDeMean, M, mr_parms);
    
    [~, M, mr_parms]  = load_mgh(fullfile(rh_De_surfs{2,1},rh_De_surfs{1,1}));
    fpRhDeMean = fullfile(dpMean, ['rh.De.' num2str(str2double(extract(lh_De_surfs{1,order(dep)}, digitsPattern))) '.fsavg.mean.mgz']);
    save_mgh(rh_De_mean, fpRhDeMean, M, mr_parms);
    
    [~, M, mr_parms]  = load_mgh(fullfile(lh_Din_surfs{2,1},lh_Din_surfs{1,1}));
    fpLhDinMean = fullfile(dpMean, ['lh.Din.' num2str(str2double(extract(lh_De_surfs{1,order(dep)}, digitsPattern))) '.fsavg.mean.mgz']);
    save_mgh(lh_Din_mean, fpLhDinMean, M, mr_parms);
    
    [~, M, mr_parms]  = load_mgh(fullfile(rh_Din_surfs{2,1},rh_Din_surfs{1,1}));
    fpRhDinMean = fullfile(dpMean, ['rh.Din.' num2str(str2double(extract(lh_De_surfs{1,order(dep)}, digitsPattern))) '.fsavg.mean.mgz']);
    save_mgh(rh_Din_mean, fpRhDinMean, M, mr_parms);
    
    [~, M, mr_parms]  = load_mgh(fullfile(lh_fec_surfs{2,1},lh_fec_surfs{1,1}));
    fpLhfecMean = fullfile(dpMean, ['lh.fextra.' num2str(str2double(extract(lh_De_surfs{1,order(dep)}, digitsPattern))) '.fsavg.mean.mgz']);
    save_mgh(lh_fec_mean, fpLhfecMean, M, mr_parms);
    
    [~, M, mr_parms]  = load_mgh(fullfile(rh_fec_surfs{2,1},rh_fec_surfs{1,1}));
    fpRhfecMean = fullfile(dpMean, ['rh.fextra.' num2str(str2double(extract(lh_De_surfs{1,order(dep)}, digitsPattern))) '.fsavg.mean.mgz']);
    save_mgh(rh_fec_mean, fpRhfecMean, M, mr_parms);
    
    [~, M, mr_parms]  = load_mgh(fullfile(lh_fin_surfs{2,1},lh_fin_surfs{1,1}));
    fpLhfinMean = fullfile(dpMean, ['lh.fneurite.' num2str(str2double(extract(lh_De_surfs{1,order(dep)}, digitsPattern))) '.fsavg.mean.mgz']);
    save_mgh(lh_fin_mean, fpLhfinMean, M, mr_parms);
    
    [~, M, mr_parms]  = load_mgh(fullfile(rh_fin_surfs{2,1},rh_fin_surfs{1,1}));
    fpRhfinMean = fullfile(dpMean, ['rh.fneurite.' num2str(str2double(extract(lh_De_surfs{1,order(dep)}, digitsPattern))) '.fsavg.mean.mgz']);
    save_mgh(rh_fin_mean, fpRhfinMean, M, mr_parms);
    
    [~, M, mr_parms]  = load_mgh(fullfile(lh_fis_surfs{2,1},lh_fis_surfs{1,1}));
    fpLhfisMean = fullfile(dpMean, ['lh.fsoma.' num2str(str2double(extract(lh_De_surfs{1,order(dep)}, digitsPattern))) '.fsavg.mean.mgz']);
    save_mgh(lh_fis_mean, fpLhfisMean, M, mr_parms);
    
    [~, M, mr_parms]  = load_mgh(fullfile(rh_fis_surfs{2,1},rh_fis_surfs{1,1}));
    fpRhfisMean = fullfile(dpMean, ['rh.fsoma.' num2str(str2double(extract(lh_De_surfs{1,order(dep)}, digitsPattern))) '.fsavg.mean.mgz']);
    save_mgh(rh_fis_mean, fpRhfisMean, M, mr_parms);
    
    [~, M, mr_parms]  = load_mgh(fullfile(lh_Rs_surfs{2,1},lh_Rs_surfs{1,1}));
    fpLhRsMean = fullfile(dpMean, ['lh.Rsoma.' num2str(str2double(extract(lh_De_surfs{1,order(dep)}, digitsPattern))) '.fsavg.mean.mgz']);
    save_mgh(lh_Rs_mean, fpLhRsMean, M, mr_parms);
    
    [~, M, mr_parms]  = load_mgh(fullfile(rh_Rs_surfs{2,1},rh_Rs_surfs{1,1}));
    fpRhRsMean = fullfile(dpMean, ['rh.Rsoma.' num2str(str2double(extract(lh_De_surfs{1,order(dep)}, digitsPattern))) '.fsavg.mean.mgz']);
    save_mgh(rh_Rs_mean, fpRhRsMean, M, mr_parms);

end