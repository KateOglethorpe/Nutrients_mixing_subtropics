%Kate Oglethorpe
%Apr 2025

%This M-file plots up global figures for the paper.
close all, clear

%% Load data

load N.mat
%load chl_final.mat
%load chl.mat
load N_conc.mat
load P_conc.mat
load Z.mat
load Z_offset.mat
load N_grad_iso_lon.mat
load N_grad_iso_lat.mat
load N_grad_dia.mat
load K_dia.mat
load K_iso.mat
load w.mat
load F_iso_lon.mat
load F_iso_lat.mat
load F_dia_diff.mat
load F_dia_adv.mat
load conv_F_iso.mat
load conv_F_dia_diff.mat
load conv_F_dia_adv.mat
load wml_F_iso.mat
load therm_F_iso.mat
load wml_F_dia_diff.mat
load therm_F_dia_diff.mat
load wml_F_dia_adv.mat
load therm_F_dia_adv.mat
load gamma_n.mat
load gamma_n_z.mat
load dens_grad_dia.mat
load F_dia_dens.mat

load glob_lon.mat
load glob_lat.mat
load z_1000.mat
load w_mld_z.mat
load w_mld_dens.mat
load gyre_mask.mat

dens_grid2 = [20.99:0.02:28.81];
dens_grid = [20.98:0.02:28.82];

load NA_lat_chl.mat
load NA_lon_chl.mat
load SA_lat_chl.mat
load SA_lon_chl.mat
load I_lat_chl.mat
load I_lon_chl.mat
load NP1_lat_chl.mat
load NP1_lon_chl.mat
load NP2_lat_chl.mat
load NP2_lon_chl.mat
load SP1_lat_chl.mat
load SP1_lon_chl.mat
load SP2_lat_chl.mat
load SP2_lon_chl.mat
load NA_mask.mat
load SA_mask.mat
load SP1_mask.mat
load SP2_mask.mat
load NP1_mask.mat
load NP2_mask.mat
load I_mask.mat

%% Plot spatial extent of subtropical gyres

%load chl [mg m-3]
load chl.mat
chl(chl == 99999) = NaN;
chl = flipud(chl);

chl_Feb_paper = load('Chl_Feb_Aqua_Modis.txt');
chl_Feb_paper(chl_Feb_paper == 99999) = NaN;
chl_Feb_paper = flipud(chl_Feb_paper);

chl_Feb = load('Feb_2021.txt');
chl_Feb(chl_Feb == 99999) = NaN;
chl_Feb = flipud(chl_Feb);

%define red box:
x_lower_lim = -80.5
x_upper_lim = -10.5
y_lower_lim = 5.5
y_upper_lim = 41.5

%over chl
figure
m_proj('miller');
m_pcolor(glob_lon, glob_lat, chl)
hold on
m_line(NA_lon_chl, NA_lat_chl, 'color','k','linewi',0.3)
m_hatch(NA_lon_chl, NA_lat_chl,'single',30,2,'color','k'); 
m_line(SA_lon_chl, SA_lat_chl, 'color','k','linewi',0.3)
m_hatch(SA_lon_chl, SA_lat_chl,'single',30,2,'color','k'); 
m_line(I_lon_chl, I_lat_chl, 'color','k','linewi',0.3)
m_hatch(I_lon_chl, I_lat_chl,'single',30,2,'color','k'); 
m_line(NP1_lon_chl, NP1_lat_chl, 'color','k','linewi',0.3)
m_hatch(NP1_lon_chl, NP1_lat_chl,'single',30,2,'color','k'); 
m_line(SP1_lon_chl, SP1_lat_chl, 'color','k','linewi',0.3)
m_hatch(SP1_lon_chl, SP1_lat_chl,'single',30,2,'color','k'); 
m_line(NP2_lon_chl, NP2_lat_chl, 'color','k','linewi',0.3)
m_hatch(NP2_lon_chl, NP2_lat_chl,'single',30,2,'color','k'); 
m_line(SP2_lon_chl, SP2_lat_chl, 'color','k','linewi',0.3)
m_hatch(SP2_lon_chl, SP2_lat_chl,'single',30,2,'color','k');
m_line([x_lower_lim; x_upper_lim; x_upper_lim; x_lower_lim; x_lower_lim],[y_upper_lim; y_upper_lim; y_lower_lim; y_lower_lim; y_upper_lim],'color','r','linewi',1)
m_grid('fontsize',14);
m_coast('color','k')
set(gcf,'color','w')
c = colorbar
cmocean('algae')
caxis([0 0.6])
c.Label.String = 'Chl-a concentration (mg m^-^3)'
set(gca,'Fontsize',14)
title('Chl 2021')

% Plot chl mask:
figure
m_proj('miller');
m_pcolor(glob_lon, glob_lat, double(chl<0.1))
hold on
m_line(NA_lon_chl, NA_lat_chl, 'color','k','linewi',0.3)
m_hatch(NA_lon_chl, NA_lat_chl,'single',30,2,'color','k'); 
m_line(SA_lon_chl, SA_lat_chl, 'color','k','linewi',0.3)
m_hatch(SA_lon_chl, SA_lat_chl,'single',30,2,'color','k'); 
m_line(I_lon_chl, I_lat_chl, 'color','k','linewi',0.3)
m_hatch(I_lon_chl, I_lat_chl,'single',30,2,'color','k'); 
m_line(NP1_lon_chl, NP1_lat_chl, 'color','k','linewi',0.3)
m_hatch(NP1_lon_chl, NP1_lat_chl,'single',30,2,'color','k'); 
m_line(SP1_lon_chl, SP1_lat_chl, 'color','k','linewi',0.3)
m_hatch(SP1_lon_chl, SP1_lat_chl,'single',30,2,'color','k'); 
m_line(NP2_lon_chl, NP2_lat_chl, 'color','k','linewi',0.3)
m_hatch(NP2_lon_chl, NP2_lat_chl,'single',30,2,'color','k'); 
m_line(SP2_lon_chl, SP2_lat_chl, 'color','k','linewi',0.3)
m_hatch(SP2_lon_chl, SP2_lat_chl,'single',30,2,'color','k');
%m_line([x_lower_lim; x_upper_lim; x_upper_lim; x_lower_lim; x_lower_lim],[y_upper_lim; y_upper_lim; y_lower_lim; y_lower_lim; y_upper_lim],'color','r','linewi',1)
m_grid('fontsize',14);
m_coast('color','k')
set(gcf,'color','w')
%c = colorbar
%cmocean('algae')
%caxis([0 0.6])
%c.Label.String = ''
set(gca,'Fontsize',14)
title('Yellow marks annual mean 2021 chl-a < 0.1 mg m^-^3')


%over chl
figure
m_proj('miller');
m_pcolor(glob_lon, glob_lat, chl_Feb)
hold on
m_line(NA_lon_chl, NA_lat_chl, 'color','k','linewi',0.3)
m_hatch(NA_lon_chl, NA_lat_chl,'single',30,2,'color','k'); 
m_line(SA_lon_chl, SA_lat_chl, 'color','k','linewi',0.3)
m_hatch(SA_lon_chl, SA_lat_chl,'single',30,2,'color','k'); 
m_line(I_lon_chl, I_lat_chl, 'color','k','linewi',0.3)
m_hatch(I_lon_chl, I_lat_chl,'single',30,2,'color','k'); 
m_line(NP1_lon_chl, NP1_lat_chl, 'color','k','linewi',0.3)
m_hatch(NP1_lon_chl, NP1_lat_chl,'single',30,2,'color','k'); 
m_line(SP1_lon_chl, SP1_lat_chl, 'color','k','linewi',0.3)
m_hatch(SP1_lon_chl, SP1_lat_chl,'single',30,2,'color','k'); 
m_line(NP2_lon_chl, NP2_lat_chl, 'color','k','linewi',0.3)
m_hatch(NP2_lon_chl, NP2_lat_chl,'single',30,2,'color','k'); 
m_line(SP2_lon_chl, SP2_lat_chl, 'color','k','linewi',0.3)
m_hatch(SP2_lon_chl, SP2_lat_chl,'single',30,2,'color','k');
m_line([x_lower_lim; x_upper_lim; x_upper_lim; x_lower_lim; x_lower_lim],[y_upper_lim; y_upper_lim; y_lower_lim; y_lower_lim; y_upper_lim],'color','r','linewi',1)
m_grid('fontsize',14);
m_coast('color','k')
set(gcf,'color','w')
c = colorbar
cmocean('algae')
caxis([0 0.6])
c.Label.String = 'Chl-a concentration (mg m^-^3)'
set(gca,'Fontsize',14)
title('Chl Feb 2021')

%over chl
figure
m_proj('miller');
m_pcolor(glob_lon, glob_lat, chl_Feb_paper)
hold on
m_line(NA_lon_chl, NA_lat_chl, 'color','k','linewi',0.3)
m_hatch(NA_lon_chl, NA_lat_chl,'single',30,2,'color','k'); 
m_line(SA_lon_chl, SA_lat_chl, 'color','k','linewi',0.3)
m_hatch(SA_lon_chl, SA_lat_chl,'single',30,2,'color','k'); 
m_line(I_lon_chl, I_lat_chl, 'color','k','linewi',0.3)
m_hatch(I_lon_chl, I_lat_chl,'single',30,2,'color','k'); 
m_line(NP1_lon_chl, NP1_lat_chl, 'color','k','linewi',0.3)
m_hatch(NP1_lon_chl, NP1_lat_chl,'single',30,2,'color','k'); 
m_line(SP1_lon_chl, SP1_lat_chl, 'color','k','linewi',0.3)
m_hatch(SP1_lon_chl, SP1_lat_chl,'single',30,2,'color','k'); 
m_line(NP2_lon_chl, NP2_lat_chl, 'color','k','linewi',0.3)
m_hatch(NP2_lon_chl, NP2_lat_chl,'single',30,2,'color','k'); 
m_line(SP2_lon_chl, SP2_lat_chl, 'color','k','linewi',0.3)
m_hatch(SP2_lon_chl, SP2_lat_chl,'single',30,2,'color','k');
m_line([x_lower_lim; x_upper_lim; x_upper_lim; x_lower_lim; x_lower_lim],[y_upper_lim; y_upper_lim; y_lower_lim; y_lower_lim; y_upper_lim],'color','r','linewi',1)
m_grid('fontsize',14);
m_coast('color','k')
set(gcf,'color','w')
c = colorbar
cmocean('algae')
caxis([0 0.6])
c.Label.String = 'Chl-a concentration (mg m^-^3)'
set(gca,'Fontsize',14)
title('Chl Feb 2021 (paper)')


%% Global supply maps

figure
t = tiledlayout(4,1,'TileSpacing','Compact');%,'Padding','Compact');
nexttile(1)
% figure
% subplot(2,2,1)
m_proj('miller', 'lat',[-40 40])
m_pcolor(glob_lon, glob_lat, wml_F_iso)
hold on
m_line(NA_lon_chl, NA_lat_chl, 'color','k','linewi',0.3)
m_hatch(NA_lon_chl, NA_lat_chl,'single',15,3,'color','k'); 
m_line(SA_lon_chl, SA_lat_chl, 'color','k','linewi',0.3)
m_hatch(SA_lon_chl, SA_lat_chl,'single',15,3,'color','k'); 
m_line(I_lon_chl, I_lat_chl, 'color','k','linewi',0.3)
m_hatch(I_lon_chl, I_lat_chl,'single',15,3,'color','k'); 
m_line(NP1_lon_chl, NP1_lat_chl, 'color','k','linewi',0.3)
m_hatch(NP1_lon_chl, NP1_lat_chl,'single',15,3,'color','k'); 
m_line(SP1_lon_chl, SP1_lat_chl, 'color','k','linewi',0.3)
m_hatch(SP1_lon_chl, SP1_lat_chl,'single',15,3,'color','k'); 
m_line(NP2_lon_chl, NP2_lat_chl, 'color','k','linewi',0.3)
m_hatch(NP2_lon_chl, NP2_lat_chl,'single',15,3,'color','k'); 
m_line(SP2_lon_chl, SP2_lat_chl, 'color','k','linewi',0.3)
m_hatch(SP2_lon_chl, SP2_lat_chl,'single',15,3,'color','k'); 
m_grid('fontsize',14);
m_coast('color','k')
set(gcf,'color','w')
cmocean('balance')
caxis([-0.17 0.17])
c = colorbar;
c.Label.String = 'mol N m^-^2 yr^-^1';
title('a) Supply by F_{iso} to the winter mixed layer','Units', 'normalized', 'Position', [0.5, 1.1, 1.2])
set(gca,'Fontsize',14)
nexttile(3)
m_proj('miller', 'lat',[-40 40])
m_pcolor(glob_lon, glob_lat, therm_F_iso)
hold on
m_line(NA_lon_chl, NA_lat_chl, 'color','k','linewi',0.3)
m_hatch(NA_lon_chl, NA_lat_chl,'single',15,3,'color','k'); 
m_line(SA_lon_chl, SA_lat_chl, 'color','k','linewi',0.3)
m_hatch(SA_lon_chl, SA_lat_chl,'single',15,3,'color','k'); 
m_line(I_lon_chl, I_lat_chl, 'color','k','linewi',0.3)
m_hatch(I_lon_chl, I_lat_chl,'single',15,3,'color','k'); 
m_line(NP1_lon_chl, NP1_lat_chl, 'color','k','linewi',0.3)
m_hatch(NP1_lon_chl, NP1_lat_chl,'single',15,3,'color','k'); 
m_line(SP1_lon_chl, SP1_lat_chl, 'color','k','linewi',0.3)
m_hatch(SP1_lon_chl, SP1_lat_chl,'single',15,3,'color','k'); 
m_line(NP2_lon_chl, NP2_lat_chl, 'color','k','linewi',0.3)
m_hatch(NP2_lon_chl, NP2_lat_chl,'single',15,3,'color','k'); 
m_line(SP2_lon_chl, SP2_lat_chl, 'color','k','linewi',0.3)
m_hatch(SP2_lon_chl, SP2_lat_chl,'single',15,3,'color','k'); 
m_grid('fontsize',14);
m_coast('color','k');
set(gcf,'color','w')
cmocean('balance')
caxis([-0.6 0.6])
c = colorbar;
c.Label.String = 'mol N m^-^2 yr^-^1';
title('c) Supply by F_{iso} to the thermocline')
set(gca,'Fontsize',14)
nexttile(2)
m_proj('miller', 'lat',[-40 40])
m_pcolor(glob_lon, glob_lat, (wml_F_dia_adv + wml_F_dia_diff))
hold on
m_line(NA_lon_chl, NA_lat_chl, 'color','k','linewi',0.3)
m_hatch(NA_lon_chl, NA_lat_chl,'single',15,3,'color','k'); 
m_line(SA_lon_chl, SA_lat_chl, 'color','k','linewi',0.3)
m_hatch(SA_lon_chl, SA_lat_chl,'single',15,3,'color','k'); 
m_line(I_lon_chl, I_lat_chl, 'color','k','linewi',0.3)
m_hatch(I_lon_chl, I_lat_chl,'single',15,3,'color','k'); 
m_line(NP1_lon_chl, NP1_lat_chl, 'color','k','linewi',0.3)
m_hatch(NP1_lon_chl, NP1_lat_chl,'single',15,3,'color','k'); 
m_line(SP1_lon_chl, SP1_lat_chl, 'color','k','linewi',0.3)
m_hatch(SP1_lon_chl, SP1_lat_chl,'single',15,3,'color','k'); 
m_line(NP2_lon_chl, NP2_lat_chl, 'color','k','linewi',0.3)
m_hatch(NP2_lon_chl, NP2_lat_chl,'single',15,3,'color','k'); 
m_line(SP2_lon_chl, SP2_lat_chl, 'color','k','linewi',0.3)
m_hatch(SP2_lon_chl, SP2_lat_chl,'single',15,3,'color','k'); 
m_grid('fontsize',14);
m_coast('color','k');
set(gcf,'color','w')
cmocean('balance')
c = colorbar;
caxis([-0.1 0.1])
c = colorbar;
c.Label.String = 'mol N m^-^2 yr^-^1';
title('b) Supply by F_{dia} + w*N to the winter mixed layer','Units', 'normalized', 'Position', [0.5, 1.1, 1.2])
set(gca,'Fontsize',14)
nexttile(4)
%subplot(2,2,4)
m_proj('miller', 'lat',[-40 40])
m_pcolor(glob_lon, glob_lat, (therm_F_dia_adv+therm_F_dia_diff))
hold on
m_line(NA_lon_chl, NA_lat_chl, 'color','k','linewi',0.3)
m_hatch(NA_lon_chl, NA_lat_chl,'single',15,3,'color','k'); 
m_line(SA_lon_chl, SA_lat_chl, 'color','k','linewi',0.3)
m_hatch(SA_lon_chl, SA_lat_chl,'single',15,3,'color','k'); 
m_line(I_lon_chl, I_lat_chl, 'color','k','linewi',0.3)
m_hatch(I_lon_chl, I_lat_chl,'single',15,3,'color','k'); 
m_line(NP1_lon_chl, NP1_lat_chl, 'color','k','linewi',0.3)
m_hatch(NP1_lon_chl, NP1_lat_chl,'single',15,3,'color','k'); 
m_line(SP1_lon_chl, SP1_lat_chl, 'color','k','linewi',0.3)
m_hatch(SP1_lon_chl, SP1_lat_chl,'single',15,3,'color','k'); 
m_line(NP2_lon_chl, NP2_lat_chl, 'color','k','linewi',0.3)
m_hatch(NP2_lon_chl, NP2_lat_chl,'single',15,3,'color','k'); 
m_line(SP2_lon_chl, SP2_lat_chl, 'color','k','linewi',0.3)
m_hatch(SP2_lon_chl, SP2_lat_chl,'single',15,3,'color','k'); 
m_grid('fontsize',14);
m_coast('color','k');
set(gcf,'color','w')
cmocean('balance')
c = colorbar;
caxis([-0.1 0.1])
c = colorbar;
c.Label.String = 'mol N m^-^2 yr^-^1';
title('d) Supply by F_{dia} + w*N to the thermocline')
set(gca,'Fontsize',14)
% 
% %turbulent mixing supply
% figure
% subplot(2,2,1)
% m_proj('miller', 'lat',[-40 40])
% m_pcolor(glob_lon, glob_lat, wml_F_dia_diff)
% hold on
% m_line(NA_lon_chl, NA_lat_chl, 'color','k','linewi',0.3)
% m_hatch(NA_lon_chl, NA_lat_chl,'single',15,3,'color','k'); 
% m_line(SA_lon_chl, SA_lat_chl, 'color','k','linewi',0.3)
% m_hatch(SA_lon_chl, SA_lat_chl,'single',15,3,'color','k'); 
% m_line(I_lon_chl, I_lat_chl, 'color','k','linewi',0.3)
% m_hatch(I_lon_chl, I_lat_chl,'single',15,3,'color','k'); 
% m_line(NP1_lon_chl, NP1_lat_chl, 'color','k','linewi',0.3)
% m_hatch(NP1_lon_chl, NP1_lat_chl,'single',15,3,'color','k'); 
% m_line(SP1_lon_chl, SP1_lat_chl, 'color','k','linewi',0.3)
% m_hatch(SP1_lon_chl, SP1_lat_chl,'single',15,3,'color','k'); 
% m_line(NP2_lon_chl, NP2_lat_chl, 'color','k','linewi',0.3)
% m_hatch(NP2_lon_chl, NP2_lat_chl,'single',15,3,'color','k'); 
% m_line(SP2_lon_chl, SP2_lat_chl, 'color','k','linewi',0.3)
% m_hatch(SP2_lon_chl, SP2_lat_chl,'single',15,3,'color','k'); 
% m_grid('fontsize',9);
% m_coast('color','k');
% set(gcf,'color','w')
% cmocean('balance')
% caxis([-0.6 0.6])
% c = colorbar;
% c.Label.String = 'mol N m^-^2 yr^-^1';
% set(gca,'Fontsize',9)
% title('F_d_i_a supply to winter mixed layer')
% subplot(2,2,2)
% m_proj('miller', 'lat',[-40 40])
% m_pcolor(glob_lon, glob_lat, (wml_F_dia_adv))
% hold on
% m_line(NA_lon_chl, NA_lat_chl, 'color','k','linewi',0.3)
% m_hatch(NA_lon_chl, NA_lat_chl,'single',15,3,'color','k'); 
% m_line(SA_lon_chl, SA_lat_chl, 'color','k','linewi',0.3)
% m_hatch(SA_lon_chl, SA_lat_chl,'single',15,3,'color','k'); 
% m_line(I_lon_chl, I_lat_chl, 'color','k','linewi',0.3)
% m_hatch(I_lon_chl, I_lat_chl,'single',15,3,'color','k'); 
% m_line(NP1_lon_chl, NP1_lat_chl, 'color','k','linewi',0.3)
% m_hatch(NP1_lon_chl, NP1_lat_chl,'single',15,3,'color','k'); 
% m_line(SP1_lon_chl, SP1_lat_chl, 'color','k','linewi',0.3)
% m_hatch(SP1_lon_chl, SP1_lat_chl,'single',15,3,'color','k'); 
% m_line(NP2_lon_chl, NP2_lat_chl, 'color','k','linewi',0.3)
% m_hatch(NP2_lon_chl, NP2_lat_chl,'single',15,3,'color','k'); 
% m_line(SP2_lon_chl, SP2_lat_chl, 'color','k','linewi',0.3)
% m_hatch(SP2_lon_chl, SP2_lat_chl,'single',15,3,'color','k'); 
% m_grid('fontsize',9);
% m_coast('color','k');
% set(gcf,'color','w')
% cmocean('balance')
% c = colorbar;
% caxis([-0.05 0.05])
% c = colorbar;
% c.Label.String = 'mol N m^-^2 yr^-^1';
% set(gca,'Fontsize',9)
% title('w*N supply to winter mixed layer')
% subplot(2,2,3)
% m_proj('miller', 'lat',[-40 40])
% m_pcolor(glob_lon, glob_lat, therm_F_dia_diff)
% hold on
% m_line(NA_lon_chl, NA_lat_chl, 'color','k','linewi',0.3)
% m_hatch(NA_lon_chl, NA_lat_chl,'single',15,3,'color','k'); 
% m_line(SA_lon_chl, SA_lat_chl, 'color','k','linewi',0.3)
% m_hatch(SA_lon_chl, SA_lat_chl,'single',15,3,'color','k'); 
% m_line(I_lon_chl, I_lat_chl, 'color','k','linewi',0.3)
% m_hatch(I_lon_chl, I_lat_chl,'single',15,3,'color','k'); 
% m_line(NP1_lon_chl, NP1_lat_chl, 'color','k','linewi',0.3)
% m_hatch(NP1_lon_chl, NP1_lat_chl,'single',15,3,'color','k'); 
% m_line(SP1_lon_chl, SP1_lat_chl, 'color','k','linewi',0.3)
% m_hatch(SP1_lon_chl, SP1_lat_chl,'single',15,3,'color','k'); 
% m_line(NP2_lon_chl, NP2_lat_chl, 'color','k','linewi',0.3)
% m_hatch(NP2_lon_chl, NP2_lat_chl,'single',15,3,'color','k'); 
% m_line(SP2_lon_chl, SP2_lat_chl, 'color','k','linewi',0.3)
% m_hatch(SP2_lon_chl, SP2_lat_chl,'single',15,3,'color','k'); 
% m_grid('fontsize',9);
% m_coast('color','k');
% set(gcf,'color','w')
% cmocean('balance')
% caxis([-0.6 0.6])
% c = colorbar;
% c.Label.String = 'mol N m^-^2 yr^-^1';
% set(gca,'Fontsize',9)
% title('F_d_i_a supply to thermocline')
% subplot(2,2,4)
% m_proj('miller', 'lat',[-40 40])
% m_pcolor(glob_lon, glob_lat, (therm_F_dia_adv))
% hold on
% m_line(NA_lon_chl, NA_lat_chl, 'color','k','linewi',0.3)
% m_hatch(NA_lon_chl, NA_lat_chl,'single',15,3,'color','k'); 
% m_line(SA_lon_chl, SA_lat_chl, 'color','k','linewi',0.3)
% m_hatch(SA_lon_chl, SA_lat_chl,'single',15,3,'color','k'); 
% m_line(I_lon_chl, I_lat_chl, 'color','k','linewi',0.3)
% m_hatch(I_lon_chl, I_lat_chl,'single',15,3,'color','k'); 
% m_line(NP1_lon_chl, NP1_lat_chl, 'color','k','linewi',0.3)
% m_hatch(NP1_lon_chl, NP1_lat_chl,'single',15,3,'color','k'); 
% m_line(SP1_lon_chl, SP1_lat_chl, 'color','k','linewi',0.3)
% m_hatch(SP1_lon_chl, SP1_lat_chl,'single',15,3,'color','k'); 
% m_line(NP2_lon_chl, NP2_lat_chl, 'color','k','linewi',0.3)
% m_hatch(NP2_lon_chl, NP2_lat_chl,'single',15,3,'color','k'); 
% m_line(SP2_lon_chl, SP2_lat_chl, 'color','k','linewi',0.3)
% m_hatch(SP2_lon_chl, SP2_lat_chl,'single',15,3,'color','k'); 
% m_grid('fontsize',9);
% m_coast('color','k');
% set(gcf,'color','w')
% cmocean('balance')
% c = colorbar;
% caxis([-0.05 0.05])
% c = colorbar;
% c.Label.String = 'mol N m^-^2 yr^-^1';
% set(gca,'Fontsize',9)
% title('w*N supply to thermocline')

%Explore supply rates on the maps...

%section = (therm_F_dia_diff+therm_F_dia_adv).*(I_mask)

%lon = glob_lon%[min_lon:max_lon]%glob_lon
%lat = glob_lat%[min_lat:max_lat]%glob_lat

%get maximum and minimum
%mAx = max(section(:));
%mIn= min(section(:));
%abs_section = abs(section);
%mIn_abs = min(abs_section(:));

%[mAx_lat, mAx_lon] = find(ismember(section, max(section(:))));
%[mIn_lat, mIn_lon] = find(ismember(section, min(section(:))));
%[mIn_abs_lat, mIn_abs_lon] = find(ismember(abs_section, min(abs_section(:))));

%stats(1,1) = mAx
%stats(1,2) = lat(mAx_lat)
%stats(1,3) = lon(mAx_lon)
%stats(2,1) = mIn
%stats(2,2) = lat(mIn_lat)
%stats(2,3) = lon(mIn_lon)
stats(3,1) = mIn_abs
%stats(3,2) = lat(mIn_abs_lat)
%stats(3,3) = lon(mIn_abs_lon)

%% Get global supply bar charts (mol m-2 yr-1 but NORMALISED by area)

%integrate supplies (mol m-2 yr-1 -> mol yr-1)
area = nan(length(glob_lat), length(glob_lon)); 
for lat = 1:(length(glob_lat)-1);
    for lon = 1:(length(glob_lon)-1);
        lat_dist = sw_dist([glob_lat(lat) glob_lat(lat+1)],[glob_lon(lon) glob_lon(lon)],'km')*1000;
        lon_dist = sw_dist([glob_lat(lat) glob_lat(lat)],[glob_lon(lon) glob_lon(lon+1)],'km')*1000;
        area(lat,lon) = lat_dist*lon_dist;
    end
end

wml_F_iso_area = wml_F_iso.*area;
wml_F_dia_diff_area = wml_F_dia_diff.*area;
wml_F_dia_adv_area = wml_F_dia_adv.*area;
therm_F_iso_area = therm_F_iso.*area;
therm_F_dia_diff_area = therm_F_dia_diff.*area;
therm_F_dia_adv_area = therm_F_dia_adv.*area;

%normalised supply for each gyre (i.e., sum(supply*area) / total area )

gyre_supply_sum = nan(2, 3, 5);

gyre_supply_sum(1, 1, 1) = sum(sum(NA_mask .* wml_F_iso_area, 'omitnan')) ./ sum(sum(NA_mask .* area, 'omitnan'));
gyre_supply_sum(2, 1, 1) = sum(sum(NA_mask .* therm_F_iso_area, 'omitnan')) ./ sum(sum(NA_mask .* area, 'omitnan'));
gyre_supply_sum(1, 2, 1) = sum(sum(NA_mask .* wml_F_dia_diff_area, 'omitnan')) ./ sum(sum(NA_mask .* area, 'omitnan'));
gyre_supply_sum(2, 2, 1) = sum(sum(NA_mask .* therm_F_dia_diff_area, 'omitnan')) ./ sum(sum(NA_mask .* area, 'omitnan'));
gyre_supply_sum(1, 3, 1) = sum(sum(NA_mask .* wml_F_dia_adv_area, 'omitnan')) ./ sum(sum(NA_mask .* area, 'omitnan'));
gyre_supply_sum(2, 3, 1) = sum(sum(NA_mask .* therm_F_dia_adv_area, 'omitnan')) ./ sum(sum(NA_mask .* area, 'omitnan'));

gyre_supply_sum(1, 1, 2) = sum(sum(SA_mask .* wml_F_iso_area, 'omitnan')) ./ sum(sum(SA_mask .* area, 'omitnan'));
gyre_supply_sum(2, 1, 2) = sum(sum(SA_mask .* therm_F_iso_area, 'omitnan')) ./ sum(sum(SA_mask .* area, 'omitnan'));
gyre_supply_sum(1, 2, 2) = sum(sum(SA_mask .* wml_F_dia_diff_area, 'omitnan')) ./ sum(sum(SA_mask .* area, 'omitnan'));
gyre_supply_sum(2, 2, 2) = sum(sum(SA_mask .* therm_F_dia_diff_area, 'omitnan')) ./ sum(sum(SA_mask .* area, 'omitnan'));
gyre_supply_sum(1, 3, 2) = sum(sum(SA_mask .* wml_F_dia_adv_area, 'omitnan')) ./ sum(sum(SA_mask .* area, 'omitnan'));
gyre_supply_sum(2, 3, 2) = sum(sum(SA_mask .* therm_F_dia_adv_area, 'omitnan')) ./ sum(sum(SA_mask .* area, 'omitnan'));

gyre_supply_sum(1, 1, 3) = sum(sum(I_mask .* wml_F_iso_area, 'omitnan')) ./ sum(sum(I_mask .* area, 'omitnan'));
gyre_supply_sum(2, 1, 3) = sum(sum(I_mask .* therm_F_iso_area, 'omitnan')) ./ sum(sum(I_mask .* area, 'omitnan'));
gyre_supply_sum(1, 2, 3) = sum(sum(I_mask .* wml_F_dia_diff_area, 'omitnan')) ./ sum(sum(I_mask .* area, 'omitnan'));
gyre_supply_sum(2, 2, 3) = sum(sum(I_mask .* therm_F_dia_diff_area, 'omitnan')) ./ sum(sum(I_mask .* area, 'omitnan'));
gyre_supply_sum(1, 3, 3) = sum(sum(I_mask .* wml_F_dia_adv_area, 'omitnan')) ./ sum(sum(I_mask .* area, 'omitnan'));
gyre_supply_sum(2, 3, 3) = sum(sum(I_mask .* therm_F_dia_adv_area, 'omitnan')) ./ sum(sum(I_mask .* area, 'omitnan'));

gyre_supply_sum(1, 1, 4) = sum(sum((NP1_mask + NP2_mask) .* wml_F_iso_area, 'omitnan')) ./ sum(sum((NP1_mask + NP2_mask) .* area, 'omitnan'));
gyre_supply_sum(2, 1, 4) = sum(sum((NP1_mask + NP2_mask) .* therm_F_iso_area, 'omitnan')) ./ sum(sum((NP1_mask + NP2_mask) .* area, 'omitnan'));
gyre_supply_sum(1, 2, 4) = sum(sum((NP1_mask + NP2_mask) .* wml_F_dia_diff_area, 'omitnan')) ./ sum(sum((NP1_mask + NP2_mask) .* area, 'omitnan'));
gyre_supply_sum(2, 2, 4) = sum(sum((NP1_mask + NP2_mask) .* therm_F_dia_diff_area, 'omitnan')) ./ sum(sum((NP1_mask + NP2_mask) .* area, 'omitnan'));
gyre_supply_sum(1, 3, 4) = sum(sum((NP1_mask + NP2_mask) .* wml_F_dia_adv_area, 'omitnan')) ./ sum(sum((NP1_mask + NP2_mask) .* area, 'omitnan'));
gyre_supply_sum(2, 3, 4) = sum(sum((NP1_mask + NP2_mask) .* therm_F_dia_adv_area, 'omitnan')) ./ sum(sum((NP1_mask + NP2_mask) .* area, 'omitnan'));

gyre_supply_sum(1, 1, 5) = sum(sum((SP1_mask + SP2_mask) .* wml_F_iso_area, 'omitnan')) ./ sum(sum((SP1_mask + SP2_mask) .* area, 'omitnan'));
gyre_supply_sum(2, 1, 5) = sum(sum((SP1_mask + SP2_mask) .* therm_F_iso_area, 'omitnan')) ./ sum(sum((SP1_mask + SP2_mask) .* area, 'omitnan'));
gyre_supply_sum(1, 2, 5) = sum(sum((SP1_mask + SP2_mask) .* wml_F_dia_diff_area, 'omitnan')) ./ sum(sum((SP1_mask + SP2_mask) .* area, 'omitnan'));
gyre_supply_sum(2, 2, 5) = sum(sum((SP1_mask + SP2_mask) .* therm_F_dia_diff_area, 'omitnan')) ./ sum(sum((SP1_mask + SP2_mask) .* area, 'omitnan'));
gyre_supply_sum(1, 3, 5) = sum(sum((SP1_mask + SP2_mask) .* wml_F_dia_adv_area, 'omitnan')) ./ sum(sum((SP1_mask + SP2_mask) .* area, 'omitnan'));
gyre_supply_sum(2, 3, 5) = sum(sum((SP1_mask + SP2_mask) .* therm_F_dia_adv_area, 'omitnan')) ./ sum(sum((SP1_mask + SP2_mask) .* area, 'omitnan'));

figure
subplot(2,1,1)
X = categorical({'North Atlantic','South Atlantic','Indian','North Pacific','South Pacific'})
X = reordercats(X,{'North Atlantic','South Atlantic','Indian','North Pacific','South Pacific'});
gyre_bar = [gyre_supply_sum(1,1,1) gyre_supply_sum(1,2,1) gyre_supply_sum(1,3,1);...
    gyre_supply_sum(1,1,2) gyre_supply_sum(1,2,2) gyre_supply_sum(1,3,2);...
    gyre_supply_sum(1,1,3) gyre_supply_sum(1,2,3) gyre_supply_sum(1,3,3);...
    gyre_supply_sum(1,1,4) gyre_supply_sum(1,2,4) gyre_supply_sum(1,3,4);...
    gyre_supply_sum(1,1,5) gyre_supply_sum(1,2,5) gyre_supply_sum(1,3,5)]%.*1000
barh(X,gyre_bar,'stacked')
set(gca,'XLim',[-0.04 0.22], 'Fontsize',14)
axis ij
xlabel('mol N m^-^2 yr^-^1')
title('a) Mixed layer supply')
%title('a) Winter mixed layer')
legend('F_i_s_o', 'F_d_i_a','w*N','Location','northeast')
subplot(2,1,2)
X = categorical({'North Atlantic','South Atlantic','Indian','North Pacific','South Pacific'})
X = reordercats(X,{'North Atlantic','South Atlantic','Indian','North Pacific','South Pacific'});
gyre_bar = [gyre_supply_sum(2,1,1) gyre_supply_sum(2,2,1) gyre_supply_sum(2,3,1);...
    gyre_supply_sum(2,1,2) gyre_supply_sum(2,2,2) gyre_supply_sum(2,3,2);...
    gyre_supply_sum(2,1,3) gyre_supply_sum(2,2,3) gyre_supply_sum(2,3,3);...
    gyre_supply_sum(2,1,4) gyre_supply_sum(2,2,4) gyre_supply_sum(2,3,4);...
    gyre_supply_sum(2,1,5) gyre_supply_sum(2,2,5) gyre_supply_sum(2,3,5)]%.*1000
barh(X,gyre_bar,'stacked');
axis ij
xlabel('mol N m^-^2 yr^-^1')
title('b) Thermocline supply')
%title('b) Upper thermocline')
% legend('F_i_s_o', 'F_d_i_a','w*N','Location','northeastoutside')
set(gca,'XLim',[-0.04 0.22],'Fontsize',14)

%percentages of Fiso vs Fdia+w*N WML for each gyre:
gyre_supply_sum(1,1,1)./(gyre_supply_sum(1,1,1) + gyre_supply_sum(1,2,1) + ...
    gyre_supply_sum(1,3,1)) 
gyre_supply_sum(1,1,2)./(gyre_supply_sum(1,1,2) + gyre_supply_sum(1,2,2) + ...
    gyre_supply_sum(1,3,2)) 
gyre_supply_sum(1,1,3)./(gyre_supply_sum(1,1,3) + gyre_supply_sum(1,2,3) + ...
    gyre_supply_sum(1,3,3)) 
gyre_supply_sum(1,1,4)./(gyre_supply_sum(1,1,4) + gyre_supply_sum(1,2,4) + ...
    gyre_supply_sum(1,3,4)) 
gyre_supply_sum(1,1,5)./(gyre_supply_sum(1,1,5) + gyre_supply_sum(1,2,5) + ...
    gyre_supply_sum(1,3,5)) 

%On average, eddy stirring provides about 50% of N to the WML. 

%percentages of Fiso vs Fdia+w*N THERMOCLINE for each gyre:
gyre_supply_sum(2,1,1)./(gyre_supply_sum(2,1,1) + gyre_supply_sum(2,2,1) + ...
    gyre_supply_sum(2,3,1)) 
gyre_supply_sum(2,1,2)./(gyre_supply_sum(2,1,2) + gyre_supply_sum(2,2,2) + ...
    gyre_supply_sum(2,3,2))
gyre_supply_sum(2,1,3)./(gyre_supply_sum(2,1,3) + gyre_supply_sum(2,2,3) + ...
    gyre_supply_sum(2,3,3))
gyre_supply_sum(2,1,4)./(gyre_supply_sum(2,1,4) + gyre_supply_sum(2,2,4) + ...
    gyre_supply_sum(2,3,4)) 
gyre_supply_sum(2,1,5)./(gyre_supply_sum(2,1,5) + gyre_supply_sum(2,2,5) + ...
    gyre_supply_sum(2,3,5)) 

%get values for how much more nitrate is supplied to thermocline than mixed layer
%north atlantic: 5 times more (400 %)
%eddy stirring makes up: 80% (ML) 98% (THERM)
%sum(gyre_supply_sum(2,:,1))./sum(gyre_supply_sum(1,:,1))
%south atlantic: 1.3 times (30 %)
%sum(gyre_supply_sum(1,:,2), 'omitnan')./sum(gyre_supply_sum(1,:,2), 'omitnan')
%indian: 9.5 times (850%)
%sum(gyre_supply_sum(2,:,3), 'omitnan')./sum(gyre_supply_sum(1,:,3), 'omitnan')
%north pacific: 7.2 times (620%)
%eddy stirring makes up 68% (ML) 97% (THERM)
%sum(gyre_supply_sum(2,:,4))./sum(gyre_supply_sum(1,:,4))
%south pacific: 2 times (100%)
%eddy stirring makes up: 57% (ML) and 81% (THERM)
%sum(gyre_supply_sum(2,:,5))./sum(gyre_supply_sum(1,:,5))

%% Supply bar chart but for expanded and shrunken gyre areas

%Expanded:

expanded_gyre_supply_sum(1, 1, 1) = sum(sum(expanded_NA_mask .* wml_F_iso_area, 'omitnan')) ./ sum(sum(expanded_NA_mask .* area, 'omitnan'));
expanded_gyre_supply_sum(2, 1, 1) = sum(sum(expanded_NA_mask .* therm_F_iso_area, 'omitnan')) ./ sum(sum(expanded_NA_mask .* area, 'omitnan'));
expanded_gyre_supply_sum(1, 2, 1) = sum(sum(expanded_NA_mask .* wml_F_dia_diff_area, 'omitnan')) ./ sum(sum(expanded_NA_mask .* area, 'omitnan'));
expanded_gyre_supply_sum(2, 2, 1) = sum(sum(expanded_NA_mask .* therm_F_dia_diff_area, 'omitnan')) ./ sum(sum(expanded_NA_mask .* area, 'omitnan'));
expanded_gyre_supply_sum(1, 3, 1) = sum(sum(expanded_NA_mask .* wml_F_dia_adv_area, 'omitnan')) ./ sum(sum(expanded_NA_mask .* area, 'omitnan'));
expanded_gyre_supply_sum(2, 3, 1) = sum(sum(expanded_NA_mask .* therm_F_dia_adv_area, 'omitnan')) ./ sum(sum(expanded_NA_mask .* area, 'omitnan'));

expanded_gyre_supply_sum(1, 1, 2) = sum(sum(expanded_SA_mask .* wml_F_iso_area, 'omitnan')) ./ sum(sum(expanded_SA_mask .* area, 'omitnan'));
expanded_gyre_supply_sum(2, 1, 2) = sum(sum(expanded_SA_mask .* therm_F_iso_area, 'omitnan')) ./ sum(sum(expanded_SA_mask .* area, 'omitnan'));
expanded_gyre_supply_sum(1, 2, 2) = sum(sum(expanded_SA_mask .* wml_F_dia_diff_area, 'omitnan')) ./ sum(sum(expanded_SA_mask .* area, 'omitnan'));
expanded_gyre_supply_sum(2, 2, 2) = sum(sum(expanded_SA_mask .* therm_F_dia_diff_area, 'omitnan')) ./ sum(sum(expanded_SA_mask .* area, 'omitnan'));
expanded_gyre_supply_sum(1, 3, 2) = sum(sum(expanded_SA_mask .* wml_F_dia_adv_area, 'omitnan')) ./ sum(sum(expanded_SA_mask .* area, 'omitnan'));
expanded_gyre_supply_sum(2, 3, 2) = sum(sum(expanded_SA_mask .* therm_F_dia_adv_area, 'omitnan')) ./ sum(sum(expanded_SA_mask .* area, 'omitnan'));

expanded_gyre_supply_sum(1, 1, 3) = sum(sum(expanded_I_mask .* wml_F_iso_area, 'omitnan')) ./ sum(sum(expanded_I_mask .* area, 'omitnan'));
expanded_gyre_supply_sum(2, 1, 3) = sum(sum(expanded_I_mask .* therm_F_iso_area, 'omitnan')) ./ sum(sum(expanded_I_mask .* area, 'omitnan'));
expanded_gyre_supply_sum(1, 2, 3) = sum(sum(expanded_I_mask .* wml_F_dia_diff_area, 'omitnan')) ./ sum(sum(expanded_I_mask .* area, 'omitnan'));
expanded_gyre_supply_sum(2, 2, 3) = sum(sum(expanded_I_mask .* therm_F_dia_diff_area, 'omitnan')) ./ sum(sum(expanded_I_mask .* area, 'omitnan'));
expanded_gyre_supply_sum(1, 3, 3) = sum(sum(expanded_I_mask .* wml_F_dia_adv_area, 'omitnan')) ./ sum(sum(expanded_I_mask .* area, 'omitnan'));
expanded_gyre_supply_sum(2, 3, 3) = sum(sum(expanded_I_mask .* therm_F_dia_adv_area, 'omitnan')) ./ sum(sum(expanded_I_mask .* area, 'omitnan'));

expanded_gyre_supply_sum(1, 1, 4) = sum(sum((expanded_NP1_mask + expanded_NP2_mask) .* wml_F_iso_area, 'omitnan')) ./ sum(sum((expanded_NP1_mask + expanded_NP2_mask) .* area, 'omitnan'));
expanded_gyre_supply_sum(2, 1, 4) = sum(sum((expanded_NP1_mask + expanded_NP2_mask) .* therm_F_iso_area, 'omitnan')) ./ sum(sum((expanded_NP1_mask + expanded_NP2_mask) .* area, 'omitnan'));
expanded_gyre_supply_sum(1, 2, 4) = sum(sum((expanded_NP1_mask + expanded_NP2_mask) .* wml_F_dia_diff_area, 'omitnan')) ./ sum(sum((expanded_NP1_mask + expanded_NP2_mask) .* area, 'omitnan'));
expanded_gyre_supply_sum(2, 2, 4) = sum(sum((expanded_NP1_mask + expanded_NP2_mask) .* therm_F_dia_diff_area, 'omitnan')) ./ sum(sum((expanded_NP1_mask + expanded_NP2_mask) .* area, 'omitnan'));
expanded_gyre_supply_sum(1, 3, 4) = sum(sum((expanded_NP1_mask + expanded_NP2_mask) .* wml_F_dia_adv_area, 'omitnan')) ./ sum(sum((expanded_NP1_mask + expanded_NP2_mask) .* area, 'omitnan'));
expanded_gyre_supply_sum(2, 3, 4) = sum(sum((expanded_NP1_mask + expanded_NP2_mask) .* therm_F_dia_adv_area, 'omitnan')) ./ sum(sum((expanded_NP1_mask + expanded_NP2_mask) .* area, 'omitnan'));

expanded_gyre_supply_sum(1, 1, 5) = sum(sum((expanded_SP1_mask + expanded_SP2_mask) .* wml_F_iso_area, 'omitnan')) ./ sum(sum((expanded_SP1_mask + expanded_SP2_mask) .* area, 'omitnan'));
expanded_gyre_supply_sum(2, 1, 5) = sum(sum((expanded_SP1_mask + expanded_SP2_mask) .* therm_F_iso_area, 'omitnan')) ./ sum(sum((expanded_SP1_mask + expanded_SP2_mask) .* area, 'omitnan'));
expanded_gyre_supply_sum(1, 2, 5) = sum(sum((expanded_SP1_mask + expanded_SP2_mask) .* wml_F_dia_diff_area, 'omitnan')) ./ sum(sum((expanded_SP1_mask + expanded_SP2_mask) .* area, 'omitnan'));
expanded_gyre_supply_sum(2, 2, 5) = sum(sum((expanded_SP1_mask + expanded_SP2_mask) .* therm_F_dia_diff_area, 'omitnan')) ./ sum(sum((expanded_SP1_mask + expanded_SP2_mask) .* area, 'omitnan'));
expanded_gyre_supply_sum(1, 3, 5) = sum(sum((expanded_SP1_mask + expanded_SP2_mask) .* wml_F_dia_adv_area, 'omitnan')) ./ sum(sum((expanded_SP1_mask + expanded_SP2_mask) .* area, 'omitnan'));
expanded_gyre_supply_sum(2, 3, 5) = sum(sum((expanded_SP1_mask + expanded_SP2_mask) .* therm_F_dia_adv_area, 'omitnan')) ./ sum(sum((expanded_SP1_mask + expanded_SP2_mask) .* area, 'omitnan'));

%percentages of Fiso vs Fdia+w*N WML for each gyre:
expanded_gyre_supply_sum(1,1,1)./(expanded_gyre_supply_sum(1,1,1) + expanded_gyre_supply_sum(1,2,1) + ...
    expanded_gyre_supply_sum(1,3,1)) 
expanded_gyre_supply_sum(1,1,2)./(expanded_gyre_supply_sum(1,1,2) + expanded_gyre_supply_sum(1,2,2) + ...
    expanded_gyre_supply_sum(1,3,2)) 
expanded_gyre_supply_sum(1,1,3)./(expanded_gyre_supply_sum(1,1,3) + expanded_gyre_supply_sum(1,2,3) + ...
    expanded_gyre_supply_sum(1,3,3)) 
expanded_gyre_supply_sum(1,1,4)./(expanded_gyre_supply_sum(1,1,4) + expanded_gyre_supply_sum(1,2,4) + ...
    expanded_gyre_supply_sum(1,3,4)) 
expanded_gyre_supply_sum(1,1,5)./(expanded_gyre_supply_sum(1,1,5) + expanded_gyre_supply_sum(1,2,5) + ...
    expanded_gyre_supply_sum(1,3,5)) 

%On average, eddy stirring provides about 50% of N to the WML. 

%percentages of Fiso vs Fdia+w*N THERMOCLINE for each gyre:
expanded_gyre_supply_sum(2,1,1)./(expanded_gyre_supply_sum(2,1,1) + expanded_gyre_supply_sum(2,2,1) + ...
    expanded_gyre_supply_sum(2,3,1)) 
expanded_gyre_supply_sum(2,1,2)./(expanded_gyre_supply_sum(2,1,2) + expanded_gyre_supply_sum(2,2,2) + ...
    expanded_gyre_supply_sum(2,3,2))
expanded_gyre_supply_sum(2,1,3)./(expanded_gyre_supply_sum(2,1,3) + expanded_gyre_supply_sum(2,2,3) + ...
    expanded_gyre_supply_sum(2,3,3))
expanded_gyre_supply_sum(2,1,4)./(expanded_gyre_supply_sum(2,1,4) + expanded_gyre_supply_sum(2,2,4) + ...
    expanded_gyre_supply_sum(2,3,4)) 
expanded_gyre_supply_sum(2,1,5)./(expanded_gyre_supply_sum(2,1,5) + expanded_gyre_supply_sum(2,2,5) + ...
    expanded_gyre_supply_sum(2,3,5)) 

%Shrunken:

shrunken_gyre_supply_sum(1, 1, 1) = sum(sum(shrunken_NA_mask .* wml_F_iso_area, 'omitnan')) ./ sum(sum(shrunken_NA_mask .* area, 'omitnan'));
shrunken_gyre_supply_sum(2, 1, 1) = sum(sum(shrunken_NA_mask .* therm_F_iso_area, 'omitnan')) ./ sum(sum(shrunken_NA_mask .* area, 'omitnan'));
shrunken_gyre_supply_sum(1, 2, 1) = sum(sum(shrunken_NA_mask .* wml_F_dia_diff_area, 'omitnan')) ./ sum(sum(shrunken_NA_mask .* area, 'omitnan'));
shrunken_gyre_supply_sum(2, 2, 1) = sum(sum(shrunken_NA_mask .* therm_F_dia_diff_area, 'omitnan')) ./ sum(sum(shrunken_NA_mask .* area, 'omitnan'));
shrunken_gyre_supply_sum(1, 3, 1) = sum(sum(shrunken_NA_mask .* wml_F_dia_adv_area, 'omitnan')) ./ sum(sum(shrunken_NA_mask .* area, 'omitnan'));
shrunken_gyre_supply_sum(2, 3, 1) = sum(sum(shrunken_NA_mask .* therm_F_dia_adv_area, 'omitnan')) ./ sum(sum(shrunken_NA_mask .* area, 'omitnan'));

shrunken_gyre_supply_sum(1, 1, 2) = sum(sum(shrunken_SA_mask .* wml_F_iso_area, 'omitnan')) ./ sum(sum(shrunken_SA_mask .* area, 'omitnan'));
shrunken_gyre_supply_sum(2, 1, 2) = sum(sum(shrunken_SA_mask .* therm_F_iso_area, 'omitnan')) ./ sum(sum(shrunken_SA_mask .* area, 'omitnan'));
shrunken_gyre_supply_sum(1, 2, 2) = sum(sum(shrunken_SA_mask .* wml_F_dia_diff_area, 'omitnan')) ./ sum(sum(shrunken_SA_mask .* area, 'omitnan'));
shrunken_gyre_supply_sum(2, 2, 2) = sum(sum(shrunken_SA_mask .* therm_F_dia_diff_area, 'omitnan')) ./ sum(sum(shrunken_SA_mask .* area, 'omitnan'));
shrunken_gyre_supply_sum(1, 3, 2) = sum(sum(shrunken_SA_mask .* wml_F_dia_adv_area, 'omitnan')) ./ sum(sum(shrunken_SA_mask .* area, 'omitnan'));
shrunken_gyre_supply_sum(2, 3, 2) = sum(sum(shrunken_SA_mask .* therm_F_dia_adv_area, 'omitnan')) ./ sum(sum(shrunken_SA_mask .* area, 'omitnan'));

shrunken_gyre_supply_sum(1, 1, 3) = sum(sum(shrunken_I_mask .* wml_F_iso_area, 'omitnan')) ./ sum(sum(shrunken_I_mask .* area, 'omitnan'));
shrunken_gyre_supply_sum(2, 1, 3) = sum(sum(shrunken_I_mask .* therm_F_iso_area, 'omitnan')) ./ sum(sum(shrunken_I_mask .* area, 'omitnan'));
shrunken_gyre_supply_sum(1, 2, 3) = sum(sum(shrunken_I_mask .* wml_F_dia_diff_area, 'omitnan')) ./ sum(sum(shrunken_I_mask .* area, 'omitnan'));
shrunken_gyre_supply_sum(2, 2, 3) = sum(sum(shrunken_I_mask .* therm_F_dia_diff_area, 'omitnan')) ./ sum(sum(shrunken_I_mask .* area, 'omitnan'));
shrunken_gyre_supply_sum(1, 3, 3) = sum(sum(shrunken_I_mask .* wml_F_dia_adv_area, 'omitnan')) ./ sum(sum(shrunken_I_mask .* area, 'omitnan'));
shrunken_gyre_supply_sum(2, 3, 3) = sum(sum(shrunken_I_mask .* therm_F_dia_adv_area, 'omitnan')) ./ sum(sum(shrunken_I_mask .* area, 'omitnan'));

shrunken_gyre_supply_sum(1, 1, 4) = sum(sum((shrunken_NP1_mask + shrunken_NP2_mask) .* wml_F_iso_area, 'omitnan')) ./ sum(sum((shrunken_NP1_mask + shrunken_NP2_mask) .* area, 'omitnan'));
shrunken_gyre_supply_sum(2, 1, 4) = sum(sum((shrunken_NP1_mask + shrunken_NP2_mask) .* therm_F_iso_area, 'omitnan')) ./ sum(sum((shrunken_NP1_mask + shrunken_NP2_mask) .* area, 'omitnan'));
shrunken_gyre_supply_sum(1, 2, 4) = sum(sum((shrunken_NP1_mask + shrunken_NP2_mask) .* wml_F_dia_diff_area, 'omitnan')) ./ sum(sum((shrunken_NP1_mask + shrunken_NP2_mask) .* area, 'omitnan'));
shrunken_gyre_supply_sum(2, 2, 4) = sum(sum((shrunken_NP1_mask + shrunken_NP2_mask) .* therm_F_dia_diff_area, 'omitnan')) ./ sum(sum((shrunken_NP1_mask + shrunken_NP2_mask) .* area, 'omitnan'));
shrunken_gyre_supply_sum(1, 3, 4) = sum(sum((shrunken_NP1_mask + shrunken_NP2_mask) .* wml_F_dia_adv_area, 'omitnan')) ./ sum(sum((shrunken_NP1_mask + shrunken_NP2_mask) .* area, 'omitnan'));
shrunken_gyre_supply_sum(2, 3, 4) = sum(sum((shrunken_NP1_mask + shrunken_NP2_mask) .* therm_F_dia_adv_area, 'omitnan')) ./ sum(sum((shrunken_NP1_mask + shrunken_NP2_mask) .* area, 'omitnan'));

shrunken_gyre_supply_sum(1, 1, 5) = sum(sum((shrunken_SP1_mask + shrunken_SP2_mask) .* wml_F_iso_area, 'omitnan')) ./ sum(sum((shrunken_SP1_mask + shrunken_SP2_mask) .* area, 'omitnan'));
shrunken_gyre_supply_sum(2, 1, 5) = sum(sum((shrunken_SP1_mask + shrunken_SP2_mask) .* therm_F_iso_area, 'omitnan')) ./ sum(sum((shrunken_SP1_mask + shrunken_SP2_mask) .* area, 'omitnan'));
shrunken_gyre_supply_sum(1, 2, 5) = sum(sum((shrunken_SP1_mask + shrunken_SP2_mask) .* wml_F_dia_diff_area, 'omitnan')) ./ sum(sum((shrunken_SP1_mask + shrunken_SP2_mask) .* area, 'omitnan'));
shrunken_gyre_supply_sum(2, 2, 5) = sum(sum((shrunken_SP1_mask + shrunken_SP2_mask) .* therm_F_dia_diff_area, 'omitnan')) ./ sum(sum((shrunken_SP1_mask + shrunken_SP2_mask) .* area, 'omitnan'));
shrunken_gyre_supply_sum(1, 3, 5) = sum(sum((shrunken_SP1_mask + shrunken_SP2_mask) .* wml_F_dia_adv_area, 'omitnan')) ./ sum(sum((shrunken_SP1_mask + shrunken_SP2_mask) .* area, 'omitnan'));
shrunken_gyre_supply_sum(2, 3, 5) = sum(sum((shrunken_SP1_mask + shrunken_SP2_mask) .* therm_F_dia_adv_area, 'omitnan')) ./ sum(sum((shrunken_SP1_mask + shrunken_SP2_mask) .* area, 'omitnan'));

%percentages of Fiso vs Fdia+w*N WML for each gyre:
shrunken_gyre_supply_sum(1,1,1)./(shrunken_gyre_supply_sum(1,1,1) + shrunken_gyre_supply_sum(1,2,1) + ...
    shrunken_gyre_supply_sum(1,3,1)) 
shrunken_gyre_supply_sum(1,1,2)./(shrunken_gyre_supply_sum(1,1,2) + shrunken_gyre_supply_sum(1,2,2) + ...
    shrunken_gyre_supply_sum(1,3,2)) 
shrunken_gyre_supply_sum(1,1,3)./(shrunken_gyre_supply_sum(1,1,3) + shrunken_gyre_supply_sum(1,2,3) + ...
    shrunken_gyre_supply_sum(1,3,3)) 
shrunken_gyre_supply_sum(1,1,4)./(shrunken_gyre_supply_sum(1,1,4) + shrunken_gyre_supply_sum(1,2,4) + ...
    shrunken_gyre_supply_sum(1,3,4)) 
shrunken_gyre_supply_sum(1,1,5)./(shrunken_gyre_supply_sum(1,1,5) + shrunken_gyre_supply_sum(1,2,5) + ...
    shrunken_gyre_supply_sum(1,3,5)) 

%On average, eddy stirring provides about 50% of N to the WML. 

%percentages of Fiso vs Fdia+w*N THERMOCLINE for each gyre:
shrunken_gyre_supply_sum(2,1,1)./(shrunken_gyre_supply_sum(2,1,1) + shrunken_gyre_supply_sum(2,2,1) + ...
    shrunken_gyre_supply_sum(2,3,1)) 
shrunken_gyre_supply_sum(2,1,2)./(shrunken_gyre_supply_sum(2,1,2) + shrunken_gyre_supply_sum(2,2,2) + ...
    shrunken_gyre_supply_sum(2,3,2))
shrunken_gyre_supply_sum(2,1,3)./(shrunken_gyre_supply_sum(2,1,3) + shrunken_gyre_supply_sum(2,2,3) + ...
    shrunken_gyre_supply_sum(2,3,3))
shrunken_gyre_supply_sum(2,1,4)./(shrunken_gyre_supply_sum(2,1,4) + shrunken_gyre_supply_sum(2,2,4) + ...
    shrunken_gyre_supply_sum(2,3,4)) 
shrunken_gyre_supply_sum(2,1,5)./(shrunken_gyre_supply_sum(2,1,5) + shrunken_gyre_supply_sum(2,2,5) + ...
    shrunken_gyre_supply_sum(2,3,5)) 

% Figure Setup
figure;

% Categories for x-axis
X = categorical({'North Atlantic', 'South Atlantic', 'Indian', 'North Pacific', 'South Pacific'});
X = reordercats(X, {'North Atlantic', 'South Atlantic', 'Indian', 'North Pacific', 'South Pacific'});

% Smaller gyre size (left column)
subplot(2, 3, 1);
gyre_bar = [shrunken_gyre_supply_sum(1,1,1), shrunken_gyre_supply_sum(1,2,1), shrunken_gyre_supply_sum(1,3,1); ...
            shrunken_gyre_supply_sum(1,1,2), shrunken_gyre_supply_sum(1,2,2), shrunken_gyre_supply_sum(1,3,2); ...
            shrunken_gyre_supply_sum(1,1,3), shrunken_gyre_supply_sum(1,2,3), shrunken_gyre_supply_sum(1,3,3); ...
            shrunken_gyre_supply_sum(1,1,4), shrunken_gyre_supply_sum(1,2,4), shrunken_gyre_supply_sum(1,3,4); ...
            shrunken_gyre_supply_sum(1,1,5), shrunken_gyre_supply_sum(1,2,5), shrunken_gyre_supply_sum(1,3,5)];
barh(X, gyre_bar, 'stacked');
axis ij;
xlabel('mol N m^-^2 yr^-^1');
title('a) Mixed layer supply', 'FontWeight', 'normal');
set(gca, 'XLim', [-0.04 0.25], 'FontSize', 12);

subplot(2, 3, 4);
gyre_bar = [shrunken_gyre_supply_sum(2,1,1), shrunken_gyre_supply_sum(2,2,1), shrunken_gyre_supply_sum(2,3,1); ...
            shrunken_gyre_supply_sum(2,1,2), shrunken_gyre_supply_sum(2,2,2), shrunken_gyre_supply_sum(2,3,2); ...
            shrunken_gyre_supply_sum(2,1,3), shrunken_gyre_supply_sum(2,2,3), shrunken_gyre_supply_sum(2,3,3); ...
            shrunken_gyre_supply_sum(2,1,4), shrunken_gyre_supply_sum(2,2,4), shrunken_gyre_supply_sum(2,3,4); ...
            shrunken_gyre_supply_sum(2,1,5), shrunken_gyre_supply_sum(2,2,5), shrunken_gyre_supply_sum(2,3,5)];
barh(X, gyre_bar, 'stacked');
axis ij;
xlabel('mol N m^-^2 yr^-^1');
title('d) Thermocline supply', 'FontWeight', 'normal');
set(gca, 'XLim', [-0.04 0.25], 'FontSize', 12);

% Normal gyre size (middle column)
subplot(2, 3, 2);
gyre_bar = [gyre_supply_sum(1,1,1), gyre_supply_sum(1,2,1), gyre_supply_sum(1,3,1); ...
            gyre_supply_sum(1,1,2), gyre_supply_sum(1,2,2), gyre_supply_sum(1,3,2); ...
            gyre_supply_sum(1,1,3), gyre_supply_sum(1,2,3), gyre_supply_sum(1,3,3); ...
            gyre_supply_sum(1,1,4), gyre_supply_sum(1,2,4), gyre_supply_sum(1,3,4); ...
            gyre_supply_sum(1,1,5), gyre_supply_sum(1,2,5), gyre_supply_sum(1,3,5)];
barh(X, gyre_bar, 'stacked');
axis ij;
xlabel('mol N m^-^2 yr^-^1');
title('b) Mixed layer supply', 'FontWeight', 'normal');
set(gca, 'XLim', [-0.04 0.25], 'FontSize', 12);
yticklabels({}); % Remove y-axis labels for the middle column

subplot(2, 3, 5);
gyre_bar = [gyre_supply_sum(2,1,1), gyre_supply_sum(2,2,1), gyre_supply_sum(2,3,1); ...
            gyre_supply_sum(2,1,2), gyre_supply_sum(2,2,2), gyre_supply_sum(2,3,2); ...
            gyre_supply_sum(2,1,3), gyre_supply_sum(2,2,3), gyre_supply_sum(2,3,3); ...
            gyre_supply_sum(2,1,4), gyre_supply_sum(2,2,4), gyre_supply_sum(2,3,4); ...
            gyre_supply_sum(2,1,5), gyre_supply_sum(2,2,5), gyre_supply_sum(2,3,5)];
barh(X, gyre_bar, 'stacked');
axis ij;
xlabel('mol N m^-^2 yr^-^1');
title('e) Thermocline supply', 'FontWeight', 'normal');
set(gca, 'XLim', [-0.04 0.25], 'FontSize', 12);
yticklabels({}); % Remove y-axis labels for the middle column

% Larger gyre size (right column)
subplot(2, 3, 3);
gyre_bar = [expanded_gyre_supply_sum(1,1,1), expanded_gyre_supply_sum(1,2,1), expanded_gyre_supply_sum(1,3,1); ...
            expanded_gyre_supply_sum(1,1,2), expanded_gyre_supply_sum(1,2,2), expanded_gyre_supply_sum(1,3,2); ...
            expanded_gyre_supply_sum(1,1,3), expanded_gyre_supply_sum(1,2,3), expanded_gyre_supply_sum(1,3,3); ...
            expanded_gyre_supply_sum(1,1,4), expanded_gyre_supply_sum(1,2,4), expanded_gyre_supply_sum(1,3,4); ...
            expanded_gyre_supply_sum(1,1,5), expanded_gyre_supply_sum(1,2,5), expanded_gyre_supply_sum(1,3,5)];
barh(X, gyre_bar, 'stacked');
axis ij;
xlabel('mol N m^-^2 yr^-^1');
title('c) Mixed layer supply', 'FontWeight', 'normal');
set(gca, 'XLim', [-0.04 0.25], 'FontSize', 12);
yticklabels({}); % Remove y-axis labels for the right column
legend('F_{iso}', 'F_{dia}', 'w*N', 'Location', 'northeast');

subplot(2, 3, 6);
gyre_bar = [expanded_gyre_supply_sum(2,1,1), expanded_gyre_supply_sum(2,2,1), expanded_gyre_supply_sum(2,3,1); ...
            expanded_gyre_supply_sum(2,1,2), expanded_gyre_supply_sum(2,2,2), expanded_gyre_supply_sum(2,3,2); ...
            expanded_gyre_supply_sum(2,1,3), expanded_gyre_supply_sum(2,2,3), expanded_gyre_supply_sum(2,3,3); ...
            expanded_gyre_supply_sum(2,1,4), expanded_gyre_supply_sum(2,2,4), expanded_gyre_supply_sum(2,3,4); ...
            expanded_gyre_supply_sum(2,1,5), expanded_gyre_supply_sum(2,2,5), expanded_gyre_supply_sum(2,3,5)];
barh(X, gyre_bar, 'stacked');
axis ij;
xlabel('mol N m^-^2 yr^-^1');
title('f) Thermocline supply', 'FontWeight', 'normal');
set(gca, 'XLim', [-0.04 0.25], 'FontSize', 12);
yticklabels({}); % Remove y-axis labels for the right column

% Add overarching titles for columns
annotation('textbox', [0.14, 0.98, 0.2, 0.03], 'String', 'Smaller Gyre Size', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
annotation('textbox', [0.42, 0.98, 0.2, 0.03], 'String', 'Normal Gyre Size', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
annotation('textbox', [0.7, 0.98, 0.2, 0.03], 'String', 'Larger Gyre Size', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');

%% Sections in 3 ocean basins showing depth of supply layers

%transects of N conc and density contours
load N_z.mat
load w_mld_z.mat
figure
%Atlantic
subplot(3,1,1)
latitude_range = (glob_lat >=-60) & (glob_lat <= 60);
pcolor(glob_lat(latitude_range), z_1000,squeeze(N_z(latitude_range,find(glob_lon==-30.5),:))')
hold on
k1 = line(glob_lat(latitude_range), -squeeze(w_mld_z(latitude_range,find(glob_lon==-30.5))),'Color','r','LineWidth',1.5,'LineStyle','--')
k2=contour(glob_lat, z_1000, squeeze(gamma_n_z(:, find(glob_lon == -30.5), :))', [27 27], 'LineColor', 'blue', 'LineWidth', 1.5,'LineStyle','--');
contour(glob_lat, z_1000,squeeze(gamma_n_z(:,find(glob_lon==-30.5),:))','ShowText','on')
shading flat; grid off
xlabel('Latitude °N')
ylabel('Depth (m)')
c = colorbar
cmocean('algae')
c.Label.String = 'Nitrate conc (µmol kg^-^1)'
title('a) Atlantic Ocean: 30.5°W')
%Pacific
subplot(3,1,2)
pcolor(glob_lat(latitude_range), z_1000,squeeze(N_z((latitude_range),find(glob_lon==-170.5),:))')
hold on
k1 = line(glob_lat(latitude_range), -squeeze(w_mld_z(latitude_range,find(glob_lon==-170.5))),'Color','r','LineWidth',1.5,'LineStyle','--')
k2=contour(glob_lat, z_1000, squeeze(gamma_n_z(:, find(glob_lon == -170.5), :))', [27 27], 'LineColor', 'blue', 'LineWidth', 1.5,'LineStyle','--');
contour(glob_lat, z_1000,squeeze(gamma_n_z(:,find(glob_lon==-170.5),:))','ShowText','on')
shading flat; grid off
xlabel('Latitude °N')
ylabel('Depth (m)')
c = colorbar
cmocean('algae')
c.Label.String = 'Nitrate conc (µmol kg^-^1)'
%lgd = legend(k1,'Mixed layer','location','southeast')
title('b) Pacific Ocean: 170.5°W')
%Indian
subplot(3,1,3)
pcolor(glob_lat(latitude_range), z_1000,squeeze(N_z((latitude_range),find(glob_lon==80.5),:))')
hold on
k1 = line(glob_lat(latitude_range), -squeeze(w_mld_z(latitude_range,find(glob_lon==80.5))),'Color','r','LineWidth',1.5,'LineStyle','--')
k2=contour(glob_lat, z_1000, squeeze(gamma_n_z(:, find(glob_lon == 80.5), :))', [27 27], 'LineColor', 'blue', 'LineWidth', 1.5,'LineStyle','--');
contour(glob_lat, z_1000,squeeze(gamma_n_z(:,find(glob_lon== 80.5),:))','ShowText','on')
shading flat; grid off
xlabel('Latitude °N')
ylabel('Depth (m)')
c = colorbar
cmocean('algae')
c.Label.String = 'Nitrate conc (µmol kg^-^1)'
%lgd = legend('Winter Mixed Layer', 'Thermocline', 'Location', 'southeast');
title('c) Indian Ocean: 80.5°E')

%% Diagnostics for w*N

min_lon = -80.5%-79.5; 
max_lon = -10.5%-5.5;
min_lat = 5.5%0.5;
max_lat = 41.5%49.5;

%Nitrate with profiles marked
figure
m_proj('mercator', 'lat', [min_lat max_lat],'long', [min_lon max_lon]);
%m_pcolor(basin_lon, basin_lat, basin_N_z(:,:,4))%_area)
hold on
m_line(NA_lon_chl, NA_lat_chl, 'color','k','linewi',0.3)
m_hatch(NA_lon_chl, NA_lat_chl,'single',30,5,'color','k'); 
m_grid('fontsize',6);
m_coast('color','k')
set(gcf,'color','w')
cmocean('algae')
caxis([0 2e-3])
c = colorbar;
c.Label.String = 'Nitrate (mol N m^-^3)';
set(gca,'Fontsize',7)
% Add profile markers
%h1 = m_plot(-48.5, 31.5, 'p', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r')
h2 = m_plot(-42.5, 27.5, '^', 'MarkerSize', 8, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b')
h3 = m_plot(-55.5, 23.5, 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y')
%h4 = m_plot(-75.5, 31.5, 's', 'MarkerSize', 8, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g')
% Define the legend labels
legend_labels = {'Profile 3', 'Profile 4'};
% Create the legend
legend([h2, h3], legend_labels, 'Location', 'southeast')

% Smooth w*N profiles:

depth_original = [20.98:0.02:28.82];
depth_original_shifted = [20.99:0.02:28.81];

% Remove mixed layer values from Kdia, dense_grad_dia, F_dia_dens, and w...
for i = 1:size(K_dia, 1)  % Loop over the x-dimension
    for j = 1:size(K_dia, 2)  % Loop over the y-dimension
        
        % Find the MLD at the current (x, y) position
        mld = w_mld_dens(i, j);
        
        % Loop over the depth levels (z-dimension)
        for k = 1:size(K_dia, 3)
            % If the depth is shallower than the MLD, set the value to NaN
            if depth_original(k) < mld
                K_dia(i, j, k) = NaN;
                F_dia_dens(i,j,k) = NaN;
                dens_grad_dia(i,j,k) = NaN;
            end
        end
        
    end
end
for i = 1:size(w, 1)  % Loop over the x-dimension
    for j = 1:size(w, 2)  % Loop over the y-dimension
        
        % Find the MLD at the current (x, y) position
        mld = w_mld_dens(i, j);
        
        % Loop over the depth levels (z-dimension)
        for k = 1:size(w, 3)
            % If the depth is shallower than the MLD, set the value to NaN
            if depth_original_shifted(k) < mld
                w(i, j, k) = NaN;
            end
        end 
        
    end
end

% Original depth range with high resolution
depth_original = [20.98:0.02:28.82];
depth_original_shifted = [20.99:0.02:28.81];

% Smooth profiles (moving average)
% Profile 1
profile_smooth_Kdia_1 = movmean(squeeze(K_dia(find(glob_lat==17.5), find(glob_lon==-55.5), :)),5);
profile_smooth_dens_grad_dia_1 = movmean(squeeze(dens_grad_dia(find(glob_lat==17.5), find(glob_lon==-55.5),:)), 5);
profile_smooth_F_dia_dens_1 = movmean(squeeze(F_dia_dens(find(glob_lat==17.5), find(glob_lon==-55.5), :)), 5);
profile_smooth_w_1 = movmean(squeeze(w(find(glob_lat==17.5), find(glob_lon==-55.5), :)), 5);

% Profile 2
profile_smooth_Kdia_2 = movmean(squeeze(K_dia(find(glob_lat==28.5), find(glob_lon==-55.5), :)), 5);
profile_smooth_dens_grad_dia_2 = movmean(squeeze(dens_grad_dia(find(glob_lat==28.5), find(glob_lon==-55.5), :)), 5);
profile_smooth_F_dia_dens_2 = movmean(squeeze(F_dia_dens(find(glob_lat==28.5), find(glob_lon==-55.5), :)), 5);
profile_smooth_w_2 = movmean(squeeze(w(find(glob_lat==28.5), find(glob_lon==-55.5), :)), 5);


% Plot smoothed w*N profiles
figure;

% Kdia
subplot(1,4,1)
plot(profile_smooth_Kdia_1,depth_original, '-b','LineWidth', 1.5)
ylim([24, 27])
set(gca, 'Color', [0.9 0.9 0.9])
axis ij
grid on
xlabel('m^2 s^{-1}')
title('K_{dia}')

% dy/dz
subplot(1,4,2)
plot(profile_smooth_dens_grad_dia_1,depth_original,'-b','LineWidth', 1.5)
ylim([24, 27])
set(gca, 'Color', [0.9 0.9 0.9])
axis ij
grid on
xlabel('kg m^{-4}')
title('dγ/dz* ')

% Fdia
subplot(1,4,3)
plot(profile_smooth_F_dia_dens_1,depth_original, '-b', 'LineWidth', 1.5)
ylim([24, 27])
set(gca, 'Color', [0.9 0.9 0.9])
axis ij
grid on
xlabel('kg m^2 s^{-1}')
title('F_{dia},_γ')

% w*
subplot(1,4,4)
plot(profile_smooth_w_1,depth_original_shifted,'-b', 'LineWidth', 1.5)
ylim([24, 27])
set(gca, 'Color', [0.9 0.9 0.9])
axis ij
grid on
xlabel('m s^{-1}')
title('w*')

% % w*N
% subplot(1,5,5)
% scatter(squeeze(F_dia_adv(find(glob_lat==27.5), find(glob_lon==-42.5),:)), [20.99:0.02:28.81], 10,'x', 'b', 'LineWidth', 1.5)
% set(gca, 'Color', [0.9 0.9 0.9])
% axis ij
% grid on
% xlabel('mol N m^{-2} yr^{-1}')
% title('w*N')

% Add a common title for all subplots
sgtitle('a) North Atlantic: 17°N, 56°W');

%Profile 2

% Create a figure
figure;

% Kdia
subplot(1,4,1)
plot(profile_smooth_Kdia_2 ,depth_original, '-r', 'LineWidth', 1.5)
ylim([25.5, 27])
set(gca, 'Color', [0.9 0.9 0.9])
axis ij
grid on
ylabel ('Density (kg m^{-3})')
xlabel('m^2 s^{-1}')
title('K_{dia}')

% dy/dz
subplot(1,4,2)
plot(profile_smooth_dens_grad_dia_2 ,depth_original, '-r', 'LineWidth', 1.5)
ylim([25.5, 27])
set(gca, 'Color', [0.9 0.9 0.9])
axis ij
grid on
xlabel('kg m^{-4}')
title('dγ/dz*')

% Fdia
subplot(1,4,3)
plot(profile_smooth_F_dia_dens_2 ,depth_original, '-r', 'LineWidth', 1.5)
ylim([25.5, 27])
set(gca, 'Color', [0.9 0.9 0.9])
axis ij
grid on
xlabel('kg m^2 s^{-1}')
title('F_{dia},_γ')

% w*
subplot(1,4,4)
plot(profile_smooth_w_2 ,depth_original_shifted, '-r', 'LineWidth', 1.5)
ylim([25.5, 27])
set(gca, 'Color', [0.9 0.9 0.9])
axis ij
grid on
xlabel('m s^{-1}')
title('w*')

% % w*N
% subplot(1,5,5)
% scatter(squeeze(F_dia_adv(find(glob_lat==23.5), find(glob_lon==-55.5),:)), [20.99:0.02:28.81], 10,'x', 'r', 'LineWidth', 1.5)
% set(gca, 'Color', [0.9 0.9 0.9])
% axis ij
% grid on
% xlabel('mol N m^{-2} yr^{-1}')
% title('w*N')

% Add a common title for all subplots
sgtitle('b) North Atlantic: 29°N, 56°W');

%% Missing Kdia plot

% Loop over the x and y dimensions
for i = 1:size(basin_K_dia_z, 1)  % Loop over the x-dimension
    for j = 1:size(basin_K_dia_z, 2)  % Loop over the y-dimension
        
        % Find the MLD at the current (x, y) position
        mld = -basin_w_mld_z(i, j);
        
        % Loop over the depth levels (z-dimension)
        for k = 1:size(basin_K_dia_z, 3)
            % If the depth is shallower than the MLD, set the value to NaN
            if z_1000(k) > mld
                basin_K_dia_z(i, j, k) = NaN;
                basin_w_z(i,j,k) = NaN;
            end
        end
        
    end
end

%plot
figure
m_proj('miller','lat',[-60 60],'lon',[-180 180]);
m_pcolor(glob_lon, glob_lat, mean(K_dia_z(:,:,[1:3]),3))
hold on
m_grid('fontsize',9);
m_coast('color','k')
set(gcf,'color','w')
c = colorbar
%cmocean('thermal')
caxis([0 6e-5])
c.Label.String = 'm^2 s^-^1'
set(gca,'Fontsize',9)
title('Average diapycnal diffusivity 0-500m')

%% Figure showing how much deeper the winter mixed layer is than euphotic and summer mixed layers.

%load euphotic layer data
euph_file = 'gridded_geospatial_montly_clim_360_720_ver_0_2.nc'; %euph layer depth (m)
euph = ncread(euph_file, 'euphotic_depth');
lon_nasa = ncread(euph_file, 'lon');       %euph 
lat_nasa = ncread(euph_file, 'lat');
%restructure data
euph = permute(euph, [2, 1,3]);
%convert euph units
euph = euph*1000;  
%map euph on 1 deg resolution
lon_woa = ncread('woa18_decav_t00_01.nc', 'lon');             
lat_woa = ncread('woa18_decav_t00_01.nc', 'lat');  
euph = interp2(lon_nasa,lat_nasa', ...    
     euph(:,:,2), lon_woa, lat_woa');

%load summer mixed layer
w_mld_winter_file = 'woa18_MLD_winter.nc';          %winter mixed layer depth (m)
w_mld_winter_z = ncread(w_mld_winter_file, 'M_an'); 
w_mld_summer_file = 'woa18_MLD_summer.nc';          %winter mixed layer depth (m)
w_mld_summer_z = ncread(w_mld_summer_file, 'M_an');  
%restructure data
w_mld_winter_z = permute(w_mld_winter_z, [2, 1]); %[m]
w_mld_summer_z = permute(w_mld_summer_z, [2, 1]); %[m]

%calculate winter mixed layer (max)
w_mld_max = nan(180,360);
for lat =1:length(glob_lat);
   for lon = 1:length(glob_lon);
       if w_mld_winter_z(lat,lon) >= w_mld_summer_z(lat,lon)
           w_mld_max(lat,lon) = w_mld_winter_z(lat,lon);
       else
           w_mld_max(lat,lon) = w_mld_summer_z(lat,lon);
       end
   end
end

%calculate summer mixed layer (min)
s_mld_min = nan(180,360);
for lat =1:length(glob_lat);
   for lon = 1:length(glob_lon);
       if w_mld_winter_z(lat,lon) >= w_mld_summer_z(lat,lon)
           s_mld_min(lat,lon) = w_mld_summer_z(lat,lon);
       else
           s_mld_min(lat,lon) = w_mld_winter_z(lat,lon);
       end
   end
end



% plot depth of WML base - depth of euphotic layer base
figure;
subplot(2,1,1)
m_proj('miller','lat',[-60 60],'lon',[-180 180]);
m_pcolor(glob_lon, glob_lat, w_mld_z - euph)
hold on
m_grid('fontsize',9);
m_coast('color','k')
set(gcf,'color','w')
c = colorbar
cmocean('deep')
caxis([0 500])
c.Label.String = 'Distance (m)'
set(gca,'Fontsize',9)
title('a) Winter mixed layer - euphotic layer')

subplot(2,1,2)
m_proj('miller','lat',[-60 60],'lon',[-180 180]);
m_pcolor(glob_lon, glob_lat, w_mld_z - s_mld_min)
hold on
m_grid('fontsize',9);
m_coast('color','k')
set(gcf,'color','w')
c = colorbar
cmocean('deep')
caxis([0 500])
c.Label.String = 'Distance (m)'
set(gca,'Fontsize',9)
title('b) Winter mixed layer - summer mixed layer')
