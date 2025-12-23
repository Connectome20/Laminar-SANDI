% Main script for analyzing laminar depth profiles of SANDI intra-soma signal fraction
% across von Economo atlas
%
% Reference:
% Lee H, ..., Huang S, Communications Biology, 2026
%
% Requirements:
% - FreeSurfer installed and configured
% - FreeSurfer fsaverage surface and annotation files available
% - Economo-Koskinas cytoarchitectonic parcellation in FreeSrufer .annot format
% - Subject-specific laminar SANDI surface maps resampled to fsaverage
% - FreeSurfer MATLAB read_annotation and load_mgh functions available
% - notBoxPlot function available
%
% Author: Hansol Lee (0000-0003-2112-1197)
% Athinoula A. Martinos Center for Biomedical Imaging,
% Massachusetts General Hospital / Harvard Medical School

clear all; close all; clc;

%%

C2_sub = {'sub-002','sub-003','sub-004','sub-005','sub-006','sub-007','sub-008','sub-009','sub-012','sub-014', ...
    'sub-015','sub-016','sub-017','sub-020','sub-035','sub-036','sub-037','sub-038b','sub-039','sub-040','sub-042'};

order = [1 11 2 4 5 6 7 8 9 10 12 13 14 15 16 17 18 19 20 21 3];

C2_dproot = '/autofs/cluster/connectome2/Bay8_C2/bids';

for sn = 1:21

C2_SANDI_dir = fullfile(C2_dproot, 'derivatives/SANDI/RF_Dis2/supres/surf/surf2Fsavg', C2_sub{sn});

lh_von = '/autofs/space/virtuoso_001/users/hlee/Bay8_C2/freesurfer/fsaverage/label/lh.economo_koskinas_fsaverabe.annot';
rh_von = '/autofs/space/virtuoso_001/users/hlee/Bay8_C2/freesurfer/fsaverage/label/rh.economo_koskinas_fsaverabe.annot';

[lh_vertices, lh_label, lh_ct] = read_annotation(lh_von);
[rh_vertices, rh_label, rh_ct] = read_annotation(rh_von);

nClass = 5;
C2_class_l = cell(1,nClass);
C2_class_r = cell(1,nClass);
for k = 1:nClass
    color_int = lh_ct.table(k+1,5);
    C2_class_l{k} = find(lh_label == color_int);
    color_int = rh_ct.table(k+1,5);
    C2_class_r{k} = find(rh_label == color_int);
end

C2_fis_maps = struct2cell(dir(fullfile(C2_SANDI_dir,'*fsoma*.mgz')));

%%

for dep = 1:length(order)

    C2_fis_l_n = fullfile(C2_fis_maps{2,order(dep)},C2_fis_maps{1,order(dep)});
    C2_fis_r_n = fullfile(C2_fis_maps{2,length(order)+order(dep)},C2_fis_maps{1,length(order)+order(dep)});
    
    C2_fis_l = load_mgh(C2_fis_l_n);
    C2_fis_r = load_mgh(C2_fis_r_n);
    
    for c = 1:nClass
        C2_fis_lc(sn,dep,c) = mean(nonzeros(C2_fis_l(C2_class_l{c})));
        C2_fis_rc(sn,dep,c) = mean(nonzeros(C2_fis_r(C2_class_r{c})));
    end
end

end

%%

C2_fis_class_mean = squeeze(mean(cat(4,C2_fis_lc,C2_fis_rc),4));

%%

close all

figure
for dd = 1:21

h=notBoxPlot(C2_fis_class_mean(:,dd,1),5*dd-5,'jitter',2,'style','sdline');
set(h.semPtch,'FaceColor',[0.42, 0.05, 0.15])
set(h.semPtch,'EdgeColor',[0.42, 0.05, 0.15])
set(h.sd,'Color',[0.42, 0.05, 0.15])
set(h.mu,'Color',[0 0 0])
set(h.data,'MarkerEdgeColor',[0.65 0.65 0.65])
set(h.data,'MarkerFaceColor',[0.65 0.65 0.65])
end

set(gcf,'Position',[0 0 1600 1000])
xticks(0:5:100),xlim([-5 105]),yticks([0.2 0.35 0.5]),ylim([0.2 0.52]),set(gca,'Fontsize',30,'FontWeight','bold')
xlabel('Depth (%)','Fontsize',50,'FontWeight','bold'),ylabel('C2.0 SANDI {\it f_i_s}','Fontsize',50,'FontWeight','bold')

figure
for dd = 1:21

h=notBoxPlot(C2_fis_class_mean(:,dd,2),5*dd-5,'jitter',2,'style','sdline');
set(h.semPtch,'FaceColor',[0.00, 0.23, 0.52])
set(h.semPtch,'EdgeColor',[0.00, 0.23, 0.52])
set(h.sd,'Color',[0.00, 0.23, 0.52])
set(h.mu,'Color',[0 0 0])
set(h.data,'MarkerEdgeColor',[0.65 0.65 0.65])
set(h.data,'MarkerFaceColor',[0.65 0.65 0.65])
end

set(gcf,'Position',[0 0 1600 1000])
xticks(0:5:100),xlim([-5 105]),yticks([0.2 0.35 0.5]),ylim([0.2 0.52]),set(gca,'Fontsize',30,'FontWeight','bold')
xlabel('Depth (%)','Fontsize',50,'FontWeight','bold'),ylabel('C2.0 SANDI {\it f_i_s}','Fontsize',50,'FontWeight','bold')

figure
for dd = 1:21

h=notBoxPlot(C2_fis_class_mean(:,dd,3),5*dd-5,'jitter',2,'style','sdline');
set(h.semPtch,'FaceColor',[0.17, 0.54, 0.17])
set(h.semPtch,'EdgeColor',[0.17, 0.54, 0.17])
set(h.sd,'Color',[0.17, 0.54, 0.17])
set(h.mu,'Color',[0 0 0])
set(h.data,'MarkerEdgeColor',[0.65 0.65 0.65])
set(h.data,'MarkerFaceColor',[0.65 0.65 0.65])
end

set(gcf,'Position',[0 0 1600 1000])
xticks(0:5:100),xlim([-5 105]),yticks([0.2 0.35 0.5]),ylim([0.2 0.52]),set(gca,'Fontsize',30,'FontWeight','bold')
xlabel('Depth (%)','Fontsize',50,'FontWeight','bold'),ylabel('C2.0 SANDI {\it f_i_s}','Fontsize',50,'FontWeight','bold')

figure
for dd = 1:21

h=notBoxPlot(C2_fis_class_mean(:,dd,4),5*dd-5,'jitter',2,'style','sdline');
set(h.semPtch,'FaceColor',[0.95, 0.82, 0.00])
set(h.semPtch,'EdgeColor',[0.95, 0.82, 0.00])
set(h.sd,'Color',[0.95, 0.82, 0.00])
set(h.mu,'Color',[0 0 0])
set(h.data,'MarkerEdgeColor',[0.65 0.65 0.65])
set(h.data,'MarkerFaceColor',[0.65 0.65 0.65])
end

set(gcf,'Position',[0 0 1600 1000])
xticks(0:5:100),xlim([-5 105]),yticks([0.2 0.35 0.5]),ylim([0.2 0.52]),set(gca,'Fontsize',30,'FontWeight','bold')
xlabel('Depth (%)','Fontsize',50,'FontWeight','bold'),ylabel('C2.0 SANDI {\it f_i_s}','Fontsize',50,'FontWeight','bold')

figure
for dd = 1:21

h=notBoxPlot(C2_fis_class_mean(:,dd,5),5*dd-5,'jitter',2,'style','sdline');
set(h.semPtch,'FaceColor',[0.77, 0.13, 0.15])
set(h.semPtch,'EdgeColor',[0.77, 0.13, 0.15])
set(h.sd,'Color',[0.77, 0.13, 0.15])
set(h.mu,'Color',[0 0 0])
set(h.data,'MarkerEdgeColor',[0.65 0.65 0.65])
set(h.data,'MarkerFaceColor',[0.65 0.65 0.65])
end

set(gcf,'Position',[0 0 1600 1000])
xticks(0:5:100),xlim([-5 105]),yticks([0.2 0.35 0.5]),ylim([0.2 0.52]),set(gca,'Fontsize',30,'FontWeight','bold')
xlabel('Depth (%)','Fontsize',50,'FontWeight','bold'),ylabel('C2.0 SANDI {\it f_i_s}','Fontsize',50,'FontWeight','bold')

%%

x = linspace(0,100,21);   % 50 depth samples from 0 to 100%

[max_1, idx_1] = max(mean(C2_fis_class_mean(:,:,1),1));
peak_depth_1 = x(idx_1);

fprintf('peak: %.3f at depth %.1f%%\n', max_1, peak_depth_1);

[max_2, idx_2] = max(mean(C2_fis_class_mean(:,:,2),1));
peak_depth_2 = x(idx_2);

fprintf('peak: %.3f at depth %.1f%%\n', max_2, peak_depth_2);

[max_3, idx_3] = max(mean(C2_fis_class_mean(:,:,3),1));
peak_depth_3 = x(idx_3);

fprintf('peak: %.3f at depth %.1f%%\n', max_3, peak_depth_3);

[max_4, idx_4] = max(mean(C2_fis_class_mean(:,:,4),1));
peak_depth_4 = x(idx_4);

fprintf('peak: %.3f at depth %.1f%%\n', max_4, peak_depth_4);

[max_5, idx_5] = max(mean(C2_fis_class_mean(:,:,5),1));
peak_depth_5 = x(idx_5);

fprintf('peak: %.3f at depth %.1f%%\n', max_5, peak_depth_5);
