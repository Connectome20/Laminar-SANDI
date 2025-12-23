% Main script for comparing laminar depth profiles of SANDI intra-soma signal fraction
% between primary motor and visual cortices
%
% Reference:
% Lee H, ..., Huang S, Communications Biology, 2026
%
% Requirements:
% - FreeSurfer installed and configured
% - FreeSurfer recon-all completed for each subject
% - FreeSurfer surface and label files available
% - Subject-specific laminar SANDI surface maps in FreeSurfer .mgz format
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

C2_SANDI_dir = fullfile(C2_dproot, 'derivatives/SANDI/RF_Dis2/supres/surf/', C2_sub{sn});

C2_4al_n = 'lh.BA4a_exvivo.thresh';
C2_4ar_n = 'rh.BA4a_exvivo.thresh';
C2_4pl_n = 'lh.BA4p_exvivo.thresh';
C2_4pr_n = 'rh.BA4p_exvivo.thresh';

C2_v1l_n = 'lh.V1_exvivo.thresh';
C2_v1r_n = 'rh.V1_exvivo.thresh';
C2_v2l_n = 'lh.V2_exvivo.thresh';
C2_v2r_n = 'rh.V2_exvivo.thresh';


dproot = '/autofs/cluster/connectome2/Bay8_C2/bids';
dpFs = fullfile(dproot, C2_sub{sn}, 'anat_process');
setenv('SUBJECTS_DIR', dpFs);

C2_4al = read_label('fs',C2_4al_n);
C2_4ar = read_label('fs',C2_4ar_n);
C2_4pl = read_label('fs',C2_4pl_n);
C2_4pr = read_label('fs',C2_4pr_n);

C2_v1l = read_label('fs',C2_v1l_n);
C2_v1r = read_label('fs',C2_v1r_n);
C2_v2l = read_label('fs',C2_v2l_n);
C2_v2r = read_label('fs',C2_v2r_n);

C2_fin_maps = struct2cell(dir(fullfile(C2_SANDI_dir,'*fneurite*.mgz')));
C2_fis_maps = struct2cell(dir(fullfile(C2_SANDI_dir,'*fsoma*.mgz')));

%%

for dep = 1:length(order)

C2_fin_l_n = fullfile(C2_fin_maps{2,order(dep)},C2_fin_maps{1,order(dep)});
C2_fis_l_n = fullfile(C2_fis_maps{2,order(dep)},C2_fis_maps{1,order(dep)});

C2_fin_r_n = fullfile(C2_fin_maps{2,length(order)+order(dep)},C2_fin_maps{1,length(order)+order(dep)});
C2_fis_r_n = fullfile(C2_fis_maps{2,length(order)+order(dep)},C2_fis_maps{1,length(order)+order(dep)});

C2_fin_l = load_mgh(C2_fin_l_n);
C2_fis_l = load_mgh(C2_fis_l_n);

C2_fin_r = load_mgh(C2_fin_r_n);
C2_fis_r = load_mgh(C2_fis_r_n);

C2_fin_4al(sn,dep) = median(nonzeros(C2_fin_l(C2_4al(:,1)+1)));
C2_fis_4al(sn,dep) = median(nonzeros(C2_fis_l(C2_4al(:,1)+1)));

C2_fin_4ar(sn,dep) = median(nonzeros(C2_fin_r(C2_4ar(:,1)+1)));
C2_fis_4ar(sn,dep) = median(nonzeros(C2_fis_r(C2_4ar(:,1)+1)));

C2_fin_4pl(sn,dep) = median(nonzeros(C2_fin_l(C2_4pl(:,1)+1)));
C2_fis_4pl(sn,dep) = median(nonzeros(C2_fis_l(C2_4pl(:,1)+1)));

C2_fin_4pr(sn,dep) = median(nonzeros(C2_fin_r(C2_4pr(:,1)+1)));
C2_fis_4pr(sn,dep) = median(nonzeros(C2_fis_r(C2_4pr(:,1)+1)));

C2_fin_v1l(sn,dep) = median(nonzeros(C2_fin_l(C2_v1l(:,1)+1)));
C2_fis_v1l(sn,dep) = median(nonzeros(C2_fis_l(C2_v1l(:,1)+1)));

C2_fin_v1r(sn,dep) = median(nonzeros(C2_fin_r(C2_v1r(:,1)+1)));
C2_fis_v1r(sn,dep) = median(nonzeros(C2_fis_r(C2_v1r(:,1)+1)));

C2_fin_v2l(sn,dep) = median(nonzeros(C2_fin_l(C2_v2l(:,1)+1)));
C2_fis_v2l(sn,dep) = median(nonzeros(C2_fis_l(C2_v2l(:,1)+1)));

C2_fin_v2r(sn,dep) = median(nonzeros(C2_fin_r(C2_v2r(:,1)+1)));
C2_fis_v2r(sn,dep) = median(nonzeros(C2_fis_r(C2_v2r(:,1)+1)));

end

end

%%

C2_fin_4a = mean(cat(3,C2_fin_4al,C2_fin_4ar),3);
C2_fis_4a = mean(cat(3,C2_fis_4al,C2_fis_4ar),3);

C2_fin_4p = mean(cat(3,C2_fin_4pl,C2_fin_4pr),3);
C2_fis_4p = mean(cat(3,C2_fis_4pl,C2_fis_4pr),3);

C2_fin_4 = mean(cat(3,C2_fin_4a,C2_fin_4p),3);
C2_fis_4 = mean(cat(3,C2_fis_4a,C2_fis_4p),3);

C2_fin_v1 = mean(cat(3,C2_fin_v1l,C2_fin_v1r),3);
C2_fis_v1 = mean(cat(3,C2_fis_v1l,C2_fis_v1r),3);

C2_fin_v2 = mean(cat(3,C2_fin_v2l,C2_fin_v2r),3);
C2_fis_v2 = mean(cat(3,C2_fis_v2l,C2_fis_v2r),3);

C2_fin_v = mean(cat(3,C2_fin_v1,C2_fin_v2),3);
C2_fis_v = mean(cat(3,C2_fis_v1,C2_fis_v2),3);

%%

close all

figure
for dd = 1:21

h=notBoxPlot(C2_fis_4(:,dd),5*dd-5,'jitter',2,'style','sdline');
set(h.semPtch,'FaceColor',[0 0.4470 0.7410])
set(h.semPtch,'EdgeColor',[0 0.4470 0.7410])
set(h.sd,'Color',[0 0.4470 0.7410])
set(h.mu,'Color',[0 0 0])
set(h.data,'MarkerEdgeColor',[0.65 0.65 0.65])
set(h.data,'MarkerFaceColor',[0.65 0.65 0.65])
end

set(gcf,'Position',[0 0 1600 1000])
xticks(0:5:100),xlim([-5 105]),yticks([0.3 0.4 0.5]),ylim([0.3 0.5]),set(gca,'Fontsize',30,'FontWeight','bold')
xlabel('Depth (%)','Fontsize',50,'FontWeight','bold'),ylabel('C2.0 SANDI {\it f_i_s}','Fontsize',50,'FontWeight','bold')

figure
for dd = 1:21
h=notBoxPlot(C2_fis_v(:,dd),5*dd-5,'jitter',2,'style','sdline');
set(h.semPtch,'FaceColor',[0.6350 0.0780 0.1840])
set(h.semPtch,'EdgeColor',[0.6350 0.0780 0.1840])
set(h.sd,'Color',[0.6350 0.0780 0.1840])
set(h.mu,'Color',[0 0 0])
set(h.data,'MarkerEdgeColor',[0.65 0.65 0.65])
set(h.data,'MarkerFaceColor',[0.65 0.65 0.65])

end

set(gcf,'Position',[0 0 1600 1000])
xticks(0:5:100),xlim([-5 105]),yticks([0.3 0.4 0.5]),ylim([0.3 0.5]),set(gca,'Fontsize',30,'FontWeight','bold')
xlabel('Depth (%)','Fontsize',50,'FontWeight','bold'),ylabel('C2.0 SANDI {\it f_i_s}','Fontsize',50,'FontWeight','bold')

figure; hold on;

for dd = 1:21

    h1 = notBoxPlot(C2_fis_4(:,dd), 5*dd-5 - 0.7, 'jitter', 1, 'style', 'sdline');
    set(h1.semPtch,'FaceColor',[0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410])
    set(h1.sd,'Color',[0 0.4470 0.7410])
    set(h1.mu,'Color',[0 0 0])
    set(h1.data,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerFaceColor',[0.65 0.65 0.65])
    set(h1.data,'Marker','none')

    h2 = notBoxPlot(C2_fis_v(:,dd), 5*dd-5 + 0.7, 'jitter', 1, 'style', 'sdline');
    set(h2.semPtch,'FaceColor',[0.6350 0.0780 0.1840],'EdgeColor',[0.6350 0.0780 0.1840])
    set(h2.sd,'Color',[0.6350 0.0780 0.1840])
    set(h2.mu,'Color',[0 0 0])
    set(h2.data,'MarkerEdgeColor',[0.65 0.65 0.65],'MarkerFaceColor',[0.65 0.65 0.65])
    set(h2.data,'Marker','none')

end

set(gcf,'Position',[0 0 1600 1000])
xticks(0:5:100)
xlim([-5 105])
yticks([0.3 0.4 0.5])
ylim([0.3 0.5])
set(gca,'FontSize',30,'FontWeight','bold')
xlabel('Depth (%)','FontSize',50,'FontWeight','bold')
ylabel('C2.0 SANDI {\it f_i_s}','FontSize',50,'FontWeight','bold')
legend({'Motor','','','','Visual'},'Location','northeast','FontSize',30)

%%

x = linspace(0,100,21);   % 21 depth samples from 0 to 100%

% ---- Motor cortex (C2_fis_4) ----
[max_motor, idx_motor] = max(mean(C2_fis_4,1));
peak_depth_motor = x(idx_motor);

fprintf('Motor cortex peak: %.3f at depth %.1f%%\n', max_motor, peak_depth_motor);

% ---- Visual cortex (C2_fis_v) ----
[max_visual, idx_visual] = max(mean(C2_fis_v,1));
peak_depth_visual = x(idx_visual);

fprintf('Visual cortex peak: %.3f at depth %.1f%%\n', max_visual, peak_depth_visual);
