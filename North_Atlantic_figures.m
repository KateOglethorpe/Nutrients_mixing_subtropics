%Kate Oglethorpe
%Mar 2023

%This M-file plots up North Atlantic figures for Code_Jul_2023.m
close all, clear

%% Load data

load N.mat
load chl_final.mat
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
load F_dia_dens.mat
load dens_grad_dia.mat

load glob_lon.mat
load glob_lat.mat
load z_1000.mat
load w_mld_z.mat
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

%% Plot spatial extent of subtropical gyres

%over chl
figure
m_proj('miller');
m_pcolor(glob_lon, glob_lat, chl_final)
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
m_grid('fontsize',9);
m_coast('color','k')
set(gcf,'color','w')
c = colorbar
cmocean('algae')
caxis([0 6e4])
c.Label.String = 'Chl-a concentration (mg m^-^3)'
set(gca,'Fontsize',9)

%% Extract data from basin of choice

basin = input('Enter basin: ')

switch basin
    
    case 'NA'
        
        %define basin
        min_lon = -80.5%-79.5; 
        max_lon = -10.5%-5.5;
        min_lat = 5.5%0.5;
        max_lat = 41.5%49.5;
        %gyre_mid_lon = -53.5;
      
        %extract data from basin
        basin_K_dia = extract_basin(K_dia,min_lon,max_lon,min_lat,max_lat);
        basin_K_iso = extract_basin(K_iso,min_lon,max_lon,min_lat,max_lat);
        basin_w = extract_basin(w,min_lon,max_lon,min_lat,max_lat);
        basin_N = extract_basin(N,min_lon,max_lon,min_lat,max_lat);
        basin_N_conc = extract_basin(N_conc,min_lon,max_lon,min_lat,max_lat);
        basin_Z = extract_basin(Z,min_lon,max_lon,min_lat,max_lat);
        basin_Z_offset = extract_basin(Z_offset,min_lon,max_lon,min_lat,max_lat);
        basin_gamma_n_z = extract_basin(gamma_n_z,min_lon,max_lon,min_lat,max_lat);
        basin_N_grad_dia = extract_basin(N_grad_dia,min_lon,max_lon,min_lat,max_lat);
        basin_N_grad_iso_lon = extract_basin(N_grad_iso_lon,min_lon,max_lon,min_lat,max_lat);
        basin_N_grad_iso_lat = extract_basin(N_grad_iso_lat,min_lat,max_lon,min_lat,max_lat);
        basin_F_iso_lat = extract_basin(F_iso_lat,min_lon,max_lon,min_lat,max_lat);
        basin_F_iso_lon = extract_basin(F_iso_lon,min_lon,max_lon,min_lat,max_lat);
        basin_F_dia_adv = extract_basin(F_dia_adv,min_lon,max_lon,min_lat,max_lat);
        basin_F_dia_diff = extract_basin(F_dia_diff,min_lon,max_lon,min_lat,max_lat);
        basin_conv_F_iso = extract_basin(conv_F_iso,min_lon,max_lon,min_lat,max_lat);
        basin_conv_F_dia_diff = extract_basin(conv_F_dia_diff,min_lon,max_lon,min_lat,max_lat);
        basin_conv_F_dia_adv = extract_basin(conv_F_dia_adv,min_lon,max_lon,min_lat,max_lat);
        basin_wml_F_iso = extract_2D_basin(wml_F_iso,min_lon,max_lon,min_lat,max_lat);
        basin_therm_F_iso = extract_2D_basin(therm_F_iso,min_lon,max_lon,min_lat,max_lat);
        basin_wml_F_dia_diff = extract_2D_basin(wml_F_dia_diff ,min_lon,max_lon,min_lat,max_lat);
        basin_therm_F_dia_diff  = extract_2D_basin(therm_F_dia_diff ,min_lon,max_lon,min_lat,max_lat);
        basin_wml_F_dia_adv= extract_2D_basin(wml_F_dia_adv ,min_lon,max_lon,min_lat,max_lat);
        basin_therm_F_dia_adv = extract_2D_basin(therm_F_dia_adv ,min_lon,max_lon,min_lat,max_lat);
        basin_N_grad_iso_lat = N_grad_iso_lat(find(glob_lat>=min_lat & glob_lat<=max_lat),...
            find(glob_lon>=min_lon & glob_lon<=max_lon),:);
        basin_Z = extract_basin(Z,min_lon,max_lon,min_lat,max_lat);
        basin_lon = [min_lon:max_lon];
        basin_lat = [min_lat:max_lat];
        basin_w_mld_z = w_mld_z(find(glob_lat>=min_lat & glob_lat<=max_lat),...
            find(glob_lon>=min_lon & glob_lon<=max_lon));
        basin_gyre_mask= gyre_mask(find(glob_lat>=min_lat & glob_lat<=max_lat),...
            find(glob_lon>=min_lon & glob_lon<=max_lon));
        basin_NA_mask = gyre_mask(find(glob_lat>=min_lat & glob_lat<=max_lat),...
            find(glob_lon>=min_lon & glob_lon<=max_lon))
        basin_F_dia_dens = extract_basin(F_dia_dens,min_lon,max_lon,min_lat,max_lat);
        basin_dens_grad_dia = extract_basin(dens_grad_dia,min_lon,max_lon,min_lat,max_lat);
    case 'SA'
        
        %define basin
        min_lon = -45.5;
        max_lon = 12.5;
        min_lat = -40.5;
        max_lat = 4.5;
        gyre_mid_lon = -25.5;
        
        %extract data from basin
        basin_K_dia = extract_basin(K_dia,min_lon,max_lon,min_lat,max_lat);
        basin_K_iso = extract_basin(K_iso,min_lon,max_lon,min_lat,max_lat);
         basin_w = extract_basin(w,min_lon,max_lon,min_lat,max_lat);
        basin_N = extract_basin(N,min_lon,max_lon,min_lat,max_lat);
        basin_N = N(find(glob_lat>=min_lat & glob_lat<=max_lat),...
            find(glob_lon>=min_lon & glob_lon<=max_lon),:);
        basin_N_conc = extract_basin(N_conc,min_lon,max_lon,min_lat,max_lat);
        basin_Z = extract_basin(Z,min_lon,max_lon,min_lat,max_lat);
        basin_Z_offset = extract_basin(Z_offset,min_lon,max_lon,min_lat,max_lat);
        basin_gamma_n_z = extract_basin(gamma_n_z,min_lon,max_lon,min_lat,max_lat);
        basin_N_grad_dia = extract_basin(N_grad_dia,min_lon,max_lon,min_lat,max_lat);
        basin_N_grad_iso_lon = extract_basin(N_grad_iso_lon,min_lon,max_lon,min_lat,max_lat);
        basin_N_grad_iso_lat = extract_basin(N_grad_iso_lat,min_lat,max_lon,min_lat,max_lat);
        basin_F_iso_lat = extract_basin(F_iso_lat,min_lon,max_lon,min_lat,max_lat);
        basin_F_iso_lon = extract_basin(F_iso_lon,min_lon,max_lon,min_lat,max_lat);
        basin_F_dia_adv = extract_basin(F_dia_adv,min_lon,max_lon,min_lat,max_lat);
        basin_F_dia_diff = extract_basin(F_dia_diff,min_lon,max_lon,min_lat,max_lat);
        basin_conv_F_iso = extract_basin(conv_F_iso,min_lon,max_lon,min_lat,max_lat);
        basin_conv_F_dia_diff = extract_basin(conv_F_dia_diff,min_lon,max_lon,min_lat,max_lat);
        basin_conv_F_dia_adv = extract_basin(conv_F_dia_adv,min_lon,max_lon,min_lat,max_lat);
        basin_wml_F_iso = extract_2D_basin(wml_F_iso,min_lon,max_lon,min_lat,max_lat);
        basin_therm_F_iso = extract_2D_basin(therm_F_iso,min_lon,max_lon,min_lat,max_lat);
        basin_wml_F_dia_diff = extract_2D_basin(wml_F_dia_diff ,min_lon,max_lon,min_lat,max_lat);
        basin_therm_F_dia_diff  = extract_2D_basin(therm_F_dia_diff ,min_lon,max_lon,min_lat,max_lat);
        basin_wml_F_dia_adv= extract_2D_basin(wml_F_dia_adv ,min_lon,max_lon,min_lat,max_lat);
        basin_therm_F_dia_adv = extract_2D_basin(therm_F_dia_adv ,min_lon,max_lon,min_lat,max_lat);
        basin_N_grad_iso_lat = N_grad_iso_lat(find(glob_lat>=min_lat & glob_lat<=max_lat),...
            find(glob_lon>=min_lon & glob_lon<=max_lon),:);
        basin_Z = extract_basin(Z,min_lon,max_lon,min_lat,max_lat);
        basin_lon = [min_lon:max_lon];
        basin_lat = [min_lat:max_lat];
        basin_w_mld_z = w_mld_z(find(glob_lat>=min_lat & glob_lat<=max_lat),...
        find(glob_lon>=min_lon & glob_lon<=max_lon));
        basin_gyre_mask= gyre_mask(find(glob_lat>=min_lat & glob_lat<=max_lat),...
        find(glob_lon>=min_lon & glob_lon<=max_lon));
        
    case 'I'
        
        %define basin
        min_lon = 19.5;
        max_lon = 129.5;
        min_lat = -40.5;;
        max_lat = -0.5;
        gyre_mid_lon = 75.5;
        
        %extract data from basin
        basin_K_dia = extract_basin(K_dia,min_lon,max_lon,min_lat,max_lat);
        basin_K_iso = extract_basin(K_iso,min_lon,max_lon,min_lat,max_lat);
         basin_w = extract_basin(w,min_lon,max_lon,min_lat,max_lat);
        basin_N = extract_basin(N,min_lon,max_lon,min_lat,max_lat);
        basin_N_conc = extract_basin(N_conc,min_lon,max_lon,min_lat,max_lat);
        basin_Z = extract_basin(Z,min_lon,max_lon,min_lat,max_lat);
        basin_Z_offset = extract_basin(Z_offset,min_lon,max_lon,min_lat,max_lat);
        basin_gamma_n_z = extract_basin(gamma_n_z,min_lon,max_lon,min_lat,max_lat);
        basin_N_grad_dia = extract_basin(N_grad_dia,min_lon,max_lon,min_lat,max_lat);
        basin_N_grad_iso_lon = extract_basin(N_grad_iso_lon,min_lon,max_lon,min_lat,max_lat);
        basin_N_grad_iso_lat = extract_basin(N_grad_iso_lat,min_lat,max_lon,min_lat,max_lat);
        basin_F_iso_lat = extract_basin(F_iso_lat,min_lon,max_lon,min_lat,max_lat);
        basin_F_iso_lon = extract_basin(F_iso_lon,min_lon,max_lon,min_lat,max_lat);
        basin_F_dia_adv = extract_basin(F_dia_adv,min_lon,max_lon,min_lat,max_lat);
        basin_F_dia_diff = extract_basin(F_dia_diff,min_lon,max_lon,min_lat,max_lat);
        basin_conv_F_iso = extract_basin(conv_F_iso,min_lon,max_lon,min_lat,max_lat);
        basin_conv_F_dia_diff = extract_basin(conv_F_dia_diff,min_lon,max_lon,min_lat,max_lat);
        basin_conv_F_dia_adv = extract_basin(conv_F_dia_adv,min_lon,max_lon,min_lat,max_lat);
        basin_wml_F_iso = extract_2D_basin(wml_F_iso,min_lon,max_lon,min_lat,max_lat);
        basin_therm_F_iso = extract_2D_basin(therm_F_iso,min_lon,max_lon,min_lat,max_lat);
        basin_wml_F_dia_diff = extract_2D_basin(wml_F_dia_diff ,min_lon,max_lon,min_lat,max_lat);
        basin_therm_F_dia_diff  = extract_2D_basin(therm_F_dia_diff ,min_lon,max_lon,min_lat,max_lat);
        basin_wml_F_dia_adv= extract_2D_basin(wml_F_dia_adv ,min_lon,max_lon,min_lat,max_lat);
        basin_therm_F_dia_adv = extract_2D_basin(therm_F_dia_adv ,min_lon,max_lon,min_lat,max_lat);
        basin_N_grad_iso_lat = N_grad_iso_lat(find(glob_lat>=min_lat & glob_lat<=max_lat),...
            find(glob_lon>=min_lon & glob_lon<=max_lon),:);
        basin_Z = extract_basin(Z,min_lon,max_lon,min_lat,max_lat);
        basin_lon = [min_lon:max_lon];
        basin_lat = [min_lat:max_lat];
        basin_w_mld_z = w_mld_z(find(glob_lat>=min_lat & glob_lat<=max_lat),...
        find(glob_lon>=min_lon & glob_lon<=max_lon));
        basin_gyre_mask= gyre_mask(find(glob_lat>=min_lat & glob_lat<=max_lat),...
            find(glob_lon>=min_lon & glob_lon<=max_lon));
        
    case 'NP'
        
        %define basin
        min_lon = 114.5;
        max_lon = -104.5;
        min_lat = 0.5;
        max_lat = 49.5;
        gyre_mid_lon = 179.5;
        
        %shift to Pacific-centric
        glob_lon = [glob_lon([181:360]);glob_lon([1:180])]; 
        N = [N(:,[181:360],: ) N(:,[1:180],:)]; 
        N_conc = [N_conc(:,[181:360],: ) N_conc(:,[1:180],:)]; 
        K_iso = [K_iso(:,[181:360],:) K_iso(:,[1:180],:)];
        K_dia= [K_dia(:,[181:360],:) K_dia(:,[1:180],:)]; 
        w= [w(:,[181:360],:) w(:,[1:180],:)]; 
        Z = [Z(:,[181:360],:) Z(:,[1:180],:)]; 
        Z_offset = [Z_offset(:,[181:360],:) Z_offset(:,[1:180],:)]; 
        w_mld_dens =[w_mld_dens(:,[181:360],:) w_mld_dens(:,[1:180])]; 
        gamma_n_z = [gamma_n_z(:,[181:360],:) gamma_n_z(:,[1:180],:)]; 
        gamma_n = [gamma_n(:,[181:360],:) gamma_n(:,[1:180],:)]; 
        N_grad_dia = [N_grad_dia(:,[181:360],:) N_grad_dia(:,[1:180],:)]; 
        N_grad_iso_lat = [N_grad_iso_lat(:,[181:360],:) N_grad_iso_lat(:,[1:180],:)]; 
        N_grad_iso_lon = [N_grad_iso_lon(:,[181:360],:) N_grad_iso_lon(:,[1:180],:)]; 
        F_iso_lon = [F_iso_lon(:,[181:360],:) F_iso_lon(:,[1:180],:)]; 
        F_iso_lat = [F_iso_lat(:,[181:360],:) F_iso_lat(:,[1:180],:)]; 
        F_dia_diff =[F_dia_diff(:,[181:360],:) F_dia_diff(:,[1:180],:)]; 
        F_dia_adv =[F_dia_adv(:,[181:360],:) F_dia_adv(:,[1:180],:)]; 
        conv_F_iso =[conv_F_iso(:,[181:360],:) conv_F_iso(:,[1:180],:)]; 
        conv_F_dia_diff =[conv_F_dia_diff(:,[181:360],:) conv_F_dia_diff(:,[1:180],:)]; 
        conv_F_dia_adv =[conv_F_dia_adv(:,[181:360],:) conv_F_dia_adv(:,[1:180],:)]; 
        wml_F_iso= [wml_F_iso(:,[181:360],:) wml_F_iso(:,[1:180],:)]; 
        therm_F_iso =[therm_F_iso(:,[181:360],:) therm_F_iso(:,[1:180],:)]; 
        wml_F_dia_diff =[wml_F_dia_diff(:,[181:360],:) wml_F_dia_diff(:,[1:180],:)];
        therm_F_dia_diff = [therm_F_dia_diff(:,[181:360],:) therm_F_dia_diff(:,[1:180],:)]; 
        wml_F_dia_adv =[wml_F_dia_adv(:,[181:360],:) wml_F_dia_adv(:,[1:180],:)];
        therm_F_dia_adv  = [therm_F_dia_adv(:,[181:360],:) therm_F_dia_adv(:,[1:180],:)];
        w_mld_z = [w_mld_z(:,[181:360]) w_mld_z(:,[1:180])];
        gyre_mask = [gyre_mask(:,[181:360]) gyre_mask(:,[1:180])];

        %extract data from basin
        basin_K_dia = extract_P_basin(K_dia,min_lon,max_lon,min_lat,max_lat);
        basin_K_iso = extract_P_basin(K_iso,min_lon,max_lon,min_lat,max_lat);
        basin_w = extract_P_basin(w,min_lon,max_lon,min_lat,max_lat);
        basin_N = extract_P_basin(N,min_lon,max_lon,min_lat,max_lat);
        basin_N_conc = extract_P_basin(N_conc,min_lon,max_lon,min_lat,max_lat);
        basin_Z = extract_P_basin(Z,min_lon,max_lon,min_lat,max_lat);
        basin_Z_offset = extract_P_basin(Z_offset,min_lon,max_lon,min_lat,max_lat);
        basin_gamma_n_z = extract_P_basin(gamma_n_z,min_lon,max_lon,min_lat,max_lat);
        basin_gamma_n = extract_P_basin(gamma_n,min_lon,max_lon,min_lat,max_lat);
        basin_N_grad_dia = extract_P_basin(N_grad_dia,min_lon,max_lon,min_lat,max_lat);
        basin_N_grad_iso_lon = extract_P_basin(N_grad_iso_lon,min_lon,max_lon,min_lat,max_lat);
        %basin_N_grad_iso_lat = extract_P_basin(N_grad_iso_lat,min_lat,max_lon,min_lat,max_lat);
        basin_F_iso_lat = extract_P_basin(F_iso_lat,min_lon,max_lon,min_lat,max_lat);
        basin_F_iso_lon = extract_P_basin(F_iso_lon,min_lon,max_lon,min_lat,max_lat);
        basin_F_dia_adv = extract_P_basin(F_dia_adv,min_lon,max_lon,min_lat,max_lat);
        basin_F_dia_diff = extract_P_basin(F_dia_diff,min_lon,max_lon,min_lat,max_lat);
        basin_conv_F_iso = extract_P_basin(conv_F_iso,min_lon,max_lon,min_lat,max_lat);
        basin_conv_F_dia_diff = extract_P_basin(conv_F_dia_diff,min_lon,max_lon,min_lat,max_lat);
        basin_conv_F_dia_adv = extract_P_basin(conv_F_dia_adv,min_lon,max_lon,min_lat,max_lat);
        basin_wml_F_iso = extract_2DP_basin(wml_F_iso,min_lon,max_lon,min_lat,max_lat);
        basin_therm_F_iso = extract_2DP_basin(therm_F_iso,min_lon,max_lon,min_lat,max_lat);
        basin_wml_F_dia_diff = extract_2DP_basin(wml_F_dia_diff ,min_lon,max_lon,min_lat,max_lat);
        basin_therm_F_dia_diff  = extract_2DP_basin(therm_F_dia_diff ,min_lon,max_lon,min_lat,max_lat);
        basin_wml_F_dia_adv= extract_2DP_basin(wml_F_dia_adv ,min_lon,max_lon,min_lat,max_lat);
        basin_therm_F_dia_adv = extract_2DP_basin(therm_F_dia_adv ,min_lon,max_lon,min_lat,max_lat);
        basin_N_grad_iso_lat = N_grad_iso_lat(find(glob_lat>=min_lat & glob_lat<=max_lat),...
            find(glob_lon>=min_lon | glob_lon<=max_lon),:);
        basin_w_mld_z = w_mld_z(find(glob_lat>=min_lat & glob_lat<=max_lat),...
        find(glob_lon>=min_lon | glob_lon<=max_lon));
        basin_gyre_mask= gyre_mask(find(glob_lat>=min_lat & glob_lat<=max_lat),...
        find(glob_lon>=min_lon | glob_lon<=max_lon));
        
        %rename max_lon for extracting basin lon and lat
        max_lon = 360 - abs(max_lon);
        basin_lon = [min_lon:max_lon];
        basin_lat = [min_lat:max_lat];

    case 'SP'
        
        %define basin
        min_lon = 139.5;
        max_lon = -62.5;
        min_lat = -49.5;
        max_lat = -0.5;
        gyre_mid_lon = 360 - abs(-151.5)
        
        %shift to Pacific-centric
        glob_lon = [glob_lon([181:360]);glob_lon([1:180])]; 
        N = [N(:,[181:360],: ) N(:,[1:180],:)]; 
        N_conc = [N_conc(:,[181:360],: ) N_conc(:,[1:180],:)]; 
        K_iso = [K_iso(:,[181:360],:) K_iso(:,[1:180],:)];
        K_dia= [K_dia(:,[181:360],:) K_dia(:,[1:180],:)]; 
        w= [w(:,[181:360],:) w(:,[1:180],:)]; 
        Z = [Z(:,[181:360],:) Z(:,[1:180],:)]; 
        Z_offset = [Z_offset(:,[181:360],:) Z_offset(:,[1:180],:)]; 
        w_mld_dens =[w_mld_dens(:,[181:360],:) w_mld_dens(:,[1:180])]; 
        gamma_n_z = [gamma_n_z(:,[181:360],:) gamma_n_z(:,[1:180],:)]; 
        gamma_n = [gamma_n(:,[181:360],:) gamma_n(:,[1:180],:)]; 
        N_grad_dia = [N_grad_dia(:,[181:360],:) N_grad_dia(:,[1:180],:)]; 
        N_grad_iso_lat = [N_grad_iso_lat(:,[181:360],:) N_grad_iso_lat(:,[1:180],:)]; 
        N_grad_iso_lon = [N_grad_iso_lon(:,[181:360],:) N_grad_iso_lon(:,[1:180],:)]; 
        F_iso_lon = [F_iso_lon(:,[181:360],:) F_iso_lon(:,[1:180],:)]; 
        F_iso_lat = [F_iso_lat(:,[181:360],:) F_iso_lat(:,[1:180],:)]; 
        F_dia_diff =[F_dia_diff(:,[181:360],:) F_dia_diff(:,[1:180],:)]; 
        F_dia_adv =[F_dia_adv(:,[181:360],:) F_dia_adv(:,[1:180],:)]; 
        conv_F_iso =[conv_F_iso(:,[181:360],:) conv_F_iso(:,[1:180],:)]; 
        conv_F_dia_diff =[conv_F_dia_diff(:,[181:360],:) conv_F_dia_diff(:,[1:180],:)]; 
        conv_F_dia_adv =[conv_F_dia_adv(:,[181:360],:) conv_F_dia_adv(:,[1:180],:)]; 
        wml_F_iso= [wml_F_iso(:,[181:360],:) wml_F_iso(:,[1:180],:)]; 
        therm_F_iso =[therm_F_iso(:,[181:360],:) therm_F_iso(:,[1:180],:)]; 
        wml_F_dia_diff =[wml_F_dia_diff(:,[181:360],:) wml_F_dia_diff(:,[1:180],:)];
        therm_F_dia_diff = [therm_F_dia_diff(:,[181:360],:) therm_F_dia_diff(:,[1:180],:)]; 
        wml_F_dia_adv =[wml_F_dia_adv(:,[181:360],:) wml_F_dia_adv(:,[1:180],:)];
        therm_F_dia_adv  = [therm_F_dia_adv(:,[181:360],:) therm_F_dia_adv(:,[1:180],:)];
        w_mld_z = [w_mld_z(:,[181:360]) w_mld_z(:,[1:180])];
        gyre_mask = [gyre_mask(:,[181:360]) gyre_mask(:,[1:180])];
         
        %extract data from basin
        basin_K_dia = extract_P_basin(K_dia,min_lon,max_lon,min_lat,max_lat);
        basin_K_iso = extract_P_basin(K_iso,min_lon,max_lon,min_lat,max_lat);
        basin_w = extract_P_basin(w,min_lon,max_lon,min_lat,max_lat);
        basin_N = extract_P_basin(N,min_lon,max_lon,min_lat,max_lat);
        basin_N_conc = extract_P_basin(N_conc,min_lon,max_lon,min_lat,max_lat);
        basin_Z = extract_P_basin(Z,min_lon,max_lon,min_lat,max_lat);
        basin_Z_offset = extract_P_basin(Z_offset,min_lon,max_lon,min_lat,max_lat);
        basin_gamma_n_z = extract_P_basin(gamma_n_z,min_lon,max_lon,min_lat,max_lat);
        basin_gamma_n = extract_P_basin(gamma_n,min_lon,max_lon,min_lat,max_lat);
        basin_N_grad_dia = extract_P_basin(N_grad_dia,min_lon,max_lon,min_lat,max_lat);
        basin_N_grad_iso_lon = extract_P_basin(N_grad_iso_lon,min_lon,max_lon,min_lat,max_lat);
        %basin_N_grad_iso_lat = extract_P_basin(N_grad_iso_lat,min_lat,max_lon,min_lat,max_lat);
        basin_F_iso_lat = extract_P_basin(F_iso_lat,min_lon,max_lon,min_lat,max_lat);
        basin_F_iso_lon = extract_P_basin(F_iso_lon,min_lon,max_lon,min_lat,max_lat);
        basin_F_dia_adv = extract_P_basin(F_dia_adv,min_lon,max_lon,min_lat,max_lat);
        basin_F_dia_diff = extract_P_basin(F_dia_diff,min_lon,max_lon,min_lat,max_lat);
        basin_conv_F_iso = extract_P_basin(conv_F_iso,min_lon,max_lon,min_lat,max_lat);
        basin_conv_F_dia_diff = extract_P_basin(conv_F_dia_diff,min_lon,max_lon,min_lat,max_lat);
        basin_conv_F_dia_adv = extract_P_basin(conv_F_dia_adv,min_lon,max_lon,min_lat,max_lat);
        basin_wml_F_iso = extract_2DP_basin(wml_F_iso,min_lon,max_lon,min_lat,max_lat);
        basin_therm_F_iso = extract_2DP_basin(therm_F_iso,min_lon,max_lon,min_lat,max_lat);
        basin_wml_F_dia_diff = extract_2DP_basin(wml_F_dia_diff ,min_lon,max_lon,min_lat,max_lat);
        basin_therm_F_dia_diff  = extract_2DP_basin(therm_F_dia_diff ,min_lon,max_lon,min_lat,max_lat);
        basin_wml_F_dia_adv= extract_2DP_basin(wml_F_dia_adv ,min_lon,max_lon,min_lat,max_lat);
        basin_therm_F_dia_adv = extract_2DP_basin(therm_F_dia_adv ,min_lon,max_lon,min_lat,max_lat);
        basin_N_grad_iso_lat = N_grad_iso_lat(find(glob_lat>=min_lat & glob_lat<=max_lat),...
            find(glob_lon>=min_lon | glob_lon<=max_lon),:);
        basin_w_mld_z = w_mld_z(find(glob_lat>=min_lat & glob_lat<=max_lat),...
        find(glob_lon>=min_lon | glob_lon<=max_lon));
        basin_gyre_mask = gyre_mask(find(glob_lat>=min_lat & glob_lat<=max_lat),...
        find(glob_lon>=min_lon | glob_lon<=max_lon));
    
        %rename max_lon for extracting basin lon and lat
        max_lon = 360 - abs(max_lon);
        basin_lon = [min_lon:max_lon];
        basin_lat = [min_lat:max_lat];
end

%% Interpolate basin data onto depth space

basin_N_z = nan(length(basin_lat),length(basin_lon), length(z_1000));
for lat = 1:length(basin_lat);
    for lon = 1:length(basin_lon);
        p_orig = dens_grid';
        N_orig = squeeze(basin_N(lat,lon,:));
        p_new = squeeze(basin_gamma_n_z(lat,lon,:));
        if sum(isnan(p_new))<45 & sum(isnan(N_orig))<393;
            a = find(~isnan(N_orig));
            N_no_nans = N_orig(a);
            p_no_nans = p_orig(a);
            basin_N_z(lat,lon,:) = interp1(p_no_nans,N_no_nans,p_new);          
        else
            basin_N_z(lat,lon,:) = NaN;
        end
    end
end

basin_N_conc_z = nan(length(basin_lat),length(basin_lon), length(z_1000));
for lat = 1:length(basin_lat);
    for lon = 1:length(basin_lon);
        p_orig = dens_grid';
        N_orig = squeeze(basin_N_conc(lat,lon,:));
        p_new = squeeze(basin_gamma_n_z(lat,lon,:));
        if sum(isnan(p_new))<45 & sum(isnan(N_orig))<393;
            a = find(~isnan(N_orig));
            N_no_nans = N_orig(a);
            p_no_nans = p_orig(a);
            basin_N_conc_z(lat,lon,:) = interp1(p_no_nans,N_no_nans,p_new);          
        else
            basin_N_conc_z(lat,lon,:) = NaN;
        end
    end
end

basin_N_grad_iso_lon_z = nan(length(basin_lat),length(basin_lon), length(z_1000));
for lat = 1:length(basin_lat);
    for lon = 1:length(basin_lon);
        p_orig = dens_grid';
        N_orig = squeeze(basin_N_grad_iso_lon(lat,lon,:));
        p_new = squeeze(basin_gamma_n_z(lat,lon,:));
        if sum(isnan(p_new))<45 & sum(isnan(N_orig))<393;
            a = find(~isnan(N_orig));
            N_no_nans = N_orig(a);
            p_no_nans = p_orig(a);
            basin_N_grad_iso_lon_z(lat,lon,:) = interp1(p_no_nans,N_no_nans,p_new);          
        else
            basin_N_grad_iso_lon_z(lat,lon,:) = NaN;
        end
    end
end
basin_N_grad_iso_lat_z = nan(length(basin_lat),length(basin_lon), length(z_1000));
for lat = 1:length(basin_lat);
    for lon = 1:length(basin_lon);
        p_orig = dens_grid';
        N_orig = squeeze(basin_N_grad_iso_lat(lat,lon,:));
        p_new = squeeze(basin_gamma_n_z(lat,lon,:));
        if sum(isnan(p_new))<45 & sum(isnan(N_orig))<393;
            a = find(~isnan(N_orig));
            N_no_nans = N_orig(a);
            p_no_nans = p_orig(a);
            basin_N_grad_iso_lat_z(lat,lon,:) = interp1(p_no_nans,N_no_nans,p_new);          
        else
            basin_N_grad_iso_lat_z(lat,lon,:) = NaN;
        end
    end
end
basin_N_grad_dia_z = nan(length(basin_lat),length(basin_lon), length(z_1000));
for lat = 1:length(basin_lat);
    for lon = 1:length(basin_lon);
        p_orig = dens_grid2';
        N_orig = squeeze(basin_N_grad_dia(lat,lon,:));
        p_new = squeeze(basin_gamma_n_z(lat,lon,:));
        if sum(isnan(p_new))<45 & sum(isnan(N_orig))<391;
            a = find(~isnan(N_orig));
            N_no_nans = N_orig(a);
            p_no_nans = p_orig(a);
            basin_N_grad_dia_z(lat,lon,:) = interp1(p_no_nans,N_no_nans,p_new);          
        else
            basin_N_grad_dia_z(lat,lon,:) = NaN;
        end
    end
end
basin_K_iso_z = nan(length(basin_lat),length(basin_lon), length(z_1000));
for lat = 1:length(basin_lat);
    for lon = 1:length(basin_lon);
        p_orig = dens_grid';
        N_orig = squeeze(basin_K_iso(lat,lon,:));
        p_new = squeeze(basin_gamma_n_z(lat,lon,:));
        if sum(isnan(p_new))<45 & sum(isnan(N_orig))<393;
            a = find(~isnan(N_orig));
            N_no_nans = N_orig(a);
            p_no_nans = p_orig(a);
            basin_K_iso_z(lat,lon,:) = interp1(p_no_nans,N_no_nans,p_new);          
        else
            basin_K_iso_z(lat,lon,:) = NaN;
        end
    end
end
basin_K_dia_z = nan(length(basin_lat),length(basin_lon), length(z_1000));
for lat = 1:length(basin_lat);
    for lon = 1:length(basin_lon);
        p_orig = dens_grid';
        N_orig = squeeze(basin_K_dia(lat,lon,:));
        p_new = squeeze(basin_gamma_n_z(lat,lon,:));
        if sum(isnan(p_new))<45 & sum(isnan(N_orig))<392;
            a = find(~isnan(N_orig));
            N_no_nans = N_orig(a);
            p_no_nans = p_orig(a);
            basin_K_dia_z(lat,lon,:) = interp1(p_no_nans,N_no_nans,p_new);          
        else
            basin_K_dia_z(lat,lon,:) = NaN;
        end
    end
end
basin_w_z = nan(length(basin_lat),length(basin_lon), length(z_1000));
for lat = 1:length(basin_lat);
    for lon = 1:length(basin_lon);
        p_orig = dens_grid2';
        N_orig = squeeze(basin_w(lat,lon,:));
        p_new = squeeze(basin_gamma_n_z(lat,lon,:));
        if sum(isnan(p_new))<45 & sum(isnan(N_orig))<391;
            a = find(~isnan(N_orig));
            N_no_nans = N_orig(a);
            p_no_nans = p_orig(a);
            basin_w_z(lat,lon,:) = interp1(p_no_nans,N_no_nans,p_new);          
        else
            basin_w_z(lat,lon,:) = NaN;
        end
    end
end

%Fiso Fdia
basin_F_iso_lon_z = nan(length(basin_lat),length(basin_lon), length(z_1000));
for lat = 1:length(basin_lat);
    for lon = 1:length(basin_lon);
        p_orig = dens_grid';
        N_orig = squeeze(basin_F_iso_lon(lat,lon,:));
        p_new = squeeze(basin_gamma_n_z(lat,lon,:));
        if sum(isnan(p_new))<45 & sum(isnan(N_orig))<393;
            a = find(~isnan(N_orig));
            N_no_nans = N_orig(a);
            p_no_nans = p_orig(a);
            basin_F_iso_lon_z(lat,lon,:) = interp1(p_no_nans,N_no_nans,p_new);          
        else
            basin_F_iso_lon_z(lat,lon,:) = NaN;
        end
    end
end
basin_F_iso_lat_z = nan(length(basin_lat),length(basin_lon), length(z_1000));
for lat = 1:length(basin_lat);
    for lon = 1:length(basin_lon);
        p_orig = dens_grid';
        N_orig = squeeze(basin_F_iso_lat(lat,lon,:));
        p_new = squeeze(basin_gamma_n_z(lat,lon,:));
        if sum(isnan(p_new))<45 & sum(isnan(N_orig))<393;
            a = find(~isnan(N_orig));
            N_no_nans = N_orig(a);
            p_no_nans = p_orig(a);
            basin_F_iso_lat_z(lat,lon,:) = interp1(p_no_nans,N_no_nans,p_new);          
        else
            basin_F_iso_lat_z(lat,lon,:) = NaN;
        end
    end
end

basin_F_dia_diff_z = nan(length(basin_lat),length(basin_lon), length(z_1000));
for lat = 1:length(basin_lat);
    for lon = 1:length(basin_lon);
        p_orig = dens_grid2';
        N_orig = squeeze(basin_F_dia_diff(lat,lon,:));
        p_new = squeeze(basin_gamma_n_z(lat,lon,:));
        if sum(isnan(p_new))<45 & sum(isnan(N_orig))<391;
            a = find(~isnan(N_orig));
            N_no_nans = N_orig(a);
            p_no_nans = p_orig(a);
            basin_F_dia_diff_z(lat,lon,:) = interp1(p_no_nans,N_no_nans,p_new);          
        else
            basin_F_dia_diff_z(lat,lon,:) = NaN;
        end
    end
end

basin_F_dia_adv_z = nan(length(basin_lat),length(basin_lon), length(z_1000));
for lat = 1:length(basin_lat);
    for lon = 1:length(basin_lon);
        p_orig = dens_grid2';
        N_orig = squeeze(basin_F_dia_adv(lat,lon,:));
        p_new = squeeze(basin_gamma_n_z(lat,lon,:));
        if sum(isnan(p_new))<45 & sum(isnan(N_orig))<391;
            a = find(~isnan(N_orig));
            N_no_nans = N_orig(a);
            p_no_nans = p_orig(a);
            basin_F_dia_adv_z(lat,lon,:) = interp1(p_no_nans,N_no_nans,p_new);          
        else
            basin_F_dia_adv_z(lat,lon,:) = NaN;
        end
    end
end

basin_dens_grad_dia_z =nan(length(basin_lat),length(basin_lon), length(z_1000));
for lat = 1:length(basin_lat);
    for lon = 1:length(basin_lon);
        p_orig = [20.97:0.02:28.83]';
        N_orig = squeeze(basin_dens_grad_dia(lat,lon,:));
        p_new = squeeze(basin_gamma_n_z(lat,lon,:));
        if sum(isnan(p_new))<45 & sum(isnan(N_orig))<391;
            a = find(~isnan(N_orig));
            N_no_nans = N_orig(a);
            p_no_nans = p_orig(a);
            basin_dens_grad_dia_z(lat,lon,:) = interp1(p_no_nans,N_no_nans,p_new);          
        else
            basin_dens_grad_dia_z(lat,lon,:) = NaN;
        end
    end
end

basin_F_dia_dens_z =nan(length(basin_lat),length(basin_lon), length(z_1000));
for lat = 1:length(basin_lat);
    for lon = 1:length(basin_lon);
        p_orig = dens_grid2';
        N_orig = squeeze(basin_F_dia_dens(lat,lon,:));
        p_new = squeeze(basin_gamma_n_z(lat,lon,:));
        if sum(isnan(p_new))<45 & sum(isnan(N_orig))<392;
            a = find(~isnan(N_orig));
            N_no_nans = N_orig(a);
            p_no_nans = p_orig(a);
            basin_F_dia_dens_z(lat,lon,:) = interp1(p_no_nans,N_no_nans,p_new);          
        else
            basin_F_dia_dens_z(lat,lon,:) = NaN;
        end
    end
end

basin_conv_F_iso_z = nan(length(basin_lat),length(basin_lon), length(z_1000));
for lat = 1:length(basin_lat);
    for lon = 1:length(basin_lon);
        p_orig = dens_grid';
        N_orig = squeeze(basin_conv_F_iso(lat,lon,:));
        p_new = squeeze(basin_gamma_n_z(lat,lon,:));
        if sum(isnan(p_new))<45 & sum(isnan(N_orig))<392;
            a = find(~isnan(N_orig));
            N_no_nans = N_orig(a);
            p_no_nans = p_orig(a);
            basin_conv_F_iso_z(lat,lon,:) = interp1(p_no_nans,N_no_nans,p_new);          
        else
            basin_conv_F_iso_z(lat,lon,:) = NaN;
        end
    end
end

basin_conv_F_dia_diff_z = nan(length(basin_lat),length(basin_lon), length(z_1000));
for lat = 1:length(basin_lat);
    for lon = 1:length(basin_lon);
        p_orig = dens_grid';
        N_orig = squeeze(basin_conv_F_dia_diff(lat,lon,:));
        p_new = squeeze(basin_gamma_n_z(lat,lon,:));
        if sum(isnan(p_new))<45 & sum(isnan(N_orig))<392;
            a = find(~isnan(N_orig));
            N_no_nans = N_orig(a);
            p_no_nans = p_orig(a);
            basin_conv_F_dia_diff_z(lat,lon,:) = interp1(p_no_nans,N_no_nans,p_new);          
        else
            basin_conv_F_dia_diff_z(lat,lon,:) = NaN;
        end
    end
end

basin_conv_F_dia_adv_z = nan(length(basin_lat),length(basin_lon), length(z_1000));
for lat = 1:length(basin_lat);
    for lon = 1:length(basin_lon);
        p_orig = dens_grid';
        N_orig = squeeze(basin_conv_F_dia_adv(lat,lon,:));
        p_new = squeeze(basin_gamma_n_z(lat,lon,:));
        if sum(isnan(p_new))<45 & sum(isnan(N_orig))<392;
            a = find(~isnan(N_orig));
            N_no_nans = N_orig(a);
            p_no_nans = p_orig(a);
            basin_conv_F_dia_adv_z(lat,lon,:) = interp1(p_no_nans,N_no_nans,p_new);          
        else
            basin_conv_F_dia_adv_z(lat,lon,:) = NaN;
        end
    end
end


% Get N from mol N m-3 to microml N kg-1 [conc]
basin_N_conc_z = basin_N_z * 1e6 ./ (basin_gamma_n_z + 1000)

% transect
%subplot(2,1,2)
figure
pcolor(basin_lat, z_1000,squeeze(basin_N_conc_z(:,find(basin_lon==-55.5),:))')
hold on
k1 = line(basin_lat, -squeeze(basin_w_mld_z(:,find(basin_lon==-55.5))),'Color','r','LineWidth',2,'LineStyle','--')
k2 = line(basin_lat, basin_Z(:,find(basin_lon==-55.5),find(dens_grid==27)),'Color','b','LineWidth',2,'LineStyle','--')
contour(basin_lat, z_1000,squeeze(basin_gamma_n_z(:,find(basin_lon==-55.5),:))','ShowText','on')
shading flat; grid off
xlabel('Latitude °N')
ylabel('Depth')
xlim([10 35])
c = colorbar
cmocean('algae')
c.Label.String = 'mol N m^-^3'
caxis([0 0.03])
title('a) N')
lgd = legend([k1 k2],'Winter mixed layer','Thermocline','location','southeast')
% lgd.Title.String = 'Layer base depths'
%title('b) North Atlantic layer depths at 55.5°W')

%% 2D plots for fluxes / maps

figure
t = tiledlayout(1,2,'TileSpacing','Compact');

% Plot 1: Nitrate on 3D maps
nexttile(1)
m_proj('mercator', 'lat', [min_lat max_lat],'long', [min_lon max_lon]);
m_pcolor(basin_lon, basin_lat, squeeze(basin_N_z(:,:,find(z_1000 == -200))))
%m_line([-55.5, -55.5], [10, 35],'Color', 'k', 'LineWidth', 2)
m_line(NA_lon_chl, NA_lat_chl, 'color','k','linewi',0.3)
m_hatch(NA_lon_chl, NA_lat_chl,'single',30,5,'color','k'); 
m_grid('fontsize',10);
m_coast('color','k')
set(gcf,'color','w')
set(gca,'fontsize',10)
title('a) N at 200m')

% Plot 2: Nitrate (mol N m^-^3)
nexttile(2)
pcolor(basin_lat, z_1000, squeeze(basin_N_z(:, find(basin_lon == -55.5), :))')
hold on
k1 = line(basin_lat, -squeeze(basin_w_mld_z(:, find(basin_lon == -55.5))), 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
k2 = line(basin_lat, basin_Z(:, find(basin_lon == -55.5), find(dens_grid == 27)), 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--');
contour(basin_lat, z_1000, squeeze(basin_gamma_n_z(:, find(basin_lon == -55.5), :))', 'ShowText', 'on');
shading flat; grid off;
xlabel('Latitude °N');
ylabel('Depth');
ylim([-800 0])
xlim([10 35]);
cmocean('algae');
c1 = colorbar;
c1.Label.String = 'mol N m^-^3';
caxis([0 0.03]);
title('b) N at 56°W');
lgd = legend([k1 k2],'Winter mixed layer','Thermocline','location','southeast')

%% Plot 2D sections & maps for nutrient gradients. 

figure
t = tiledlayout(3,2,'TileSpacing','Compact');

% Plot 1: dN/dz
nexttile(1)
m_proj('mercator', 'lat', [min_lat max_lat],'long', [min_lon max_lon]);
m_pcolor(basin_lon, basin_lat, squeeze(basin_N_grad_dia_z(:,:,find(z_1000 == -200))))
%m_line([-55.5, -55.5], [10, 35],'Color', 'k', 'LineWidth', 2)
m_line(NA_lon_chl, NA_lat_chl, 'color','k','linewi',0.3)
m_hatch(NA_lon_chl, NA_lat_chl,'single',30,5,'color','k'); 
m_grid('fontsize',10);
m_coast('color','k')
set(gcf,'color','w')
set(gca,'fontsize',10)
caxis([-1e-4 1e-4]);
title('a) dN/dz^* at 200m')

% Plot 2: dN/dz
nexttile(2)
pcolor(basin_lat, z_1000, squeeze(basin_N_grad_dia_z(:, find(basin_lon == -55.5), :))')
hold on
k1 = line(basin_lat, -squeeze(basin_w_mld_z(:, find(basin_lon == -55.5))), 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
k2 = line(basin_lat, basin_Z(:, find(basin_lon == -55.5), find(dens_grid == 27)), 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--');
contour(basin_lat, z_1000, squeeze(basin_gamma_n_z(:, find(basin_lon == -55.5), :))', 'ShowText', 'on');
shading flat; grid off;
xlabel('Latitude °N');
ylabel('Depth');
ylim([-800 0])
xlim([10 35]);
cmocean('algae');
c1 = colorbar;
c1.Label.String = 'mol N m^-^4';
caxis([-1e-4 1e-4]);
title('b) dN/dz^* at 56°W');
lgd = legend([k1 k2],'Winter mixed layer','Thermocline','location','southeast')

% Plot 3: Lon N grad
nexttile(3)
m_proj('mercator', 'lat', [min_lat max_lat],'long', [min_lon max_lon]);
m_pcolor(basin_lon, basin_lat, squeeze(basin_N_grad_iso_lon_z(:,:,find(z_1000 == -200))))
%m_line([-55.5, -55.5], [10, 35],'Color', 'k', 'LineWidth', 2)
m_line(NA_lon_chl, NA_lat_chl, 'color','k','linewi',0.3)
m_hatch(NA_lon_chl, NA_lat_chl,'single',30,5,'color','k'); 
m_grid('fontsize',10);
m_coast('color','k')
set(gcf,'color','w')
set(gca,'fontsize',10)
caxis([-2e-8 2e-8]);
title('c) Lon ∇_i_s_o N at 200m')

% Plot 4: Lon N grad
nexttile(4)
pcolor(basin_lat, z_1000, squeeze(basin_N_grad_iso_lon_z(:, find(basin_lon == -55.5), :))')
hold on
k1 = line(basin_lat, -squeeze(basin_w_mld_z(:, find(basin_lon == -55.5))), 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
k2 = line(basin_lat, basin_Z(:, find(basin_lon == -55.5), find(dens_grid == 27)), 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--');
contour(basin_lat, z_1000, squeeze(basin_gamma_n_z(:, find(basin_lon == -55.5), :))', 'ShowText', 'on');
shading flat; grid off;
xlabel('Latitude °N');
ylabel('Depth');
ylim([-800 0])
xlim([10 35]);
cmocean('balance');
c1 = colorbar;
c1.Label.String = 'mol N m^-^4';
caxis([-2e-8 2e-8]);
title('d) Lon ∇_i_s_o N at 56°W');
lgd = legend([k1 k2],'Winter mixed layer','Thermocline','location','southeast')

% Plot 3: Lat N grad
nexttile(5)
m_proj('mercator', 'lat', [min_lat max_lat],'long', [min_lon max_lon]);
m_pcolor(basin_lon, basin_lat, squeeze(basin_N_grad_iso_lat_z(:,:,find(z_1000 == -200))))
%m_line([-55.5, -55.5], [10, 35],'Color', 'k', 'LineWidth', 2)
m_line(NA_lon_chl, NA_lat_chl, 'color','k','linewi',0.3)
m_hatch(NA_lon_chl, NA_lat_chl,'single',30,5,'color','k'); 
m_grid('fontsize',10);
m_coast('color','k')
set(gcf,'color','w')
set(gca,'fontsize',10)
caxis([-2e-8 2e-8]);
title('e) Lat ∇_i_s_o N  at 200m')

% Plot 4: Lat N grad
nexttile(6)
pcolor(basin_lat, z_1000, squeeze(basin_N_grad_iso_lat_z(:, find(basin_lon == -55.5), :))')
hold on
k1 = line(basin_lat, -squeeze(basin_w_mld_z(:, find(basin_lon == -55.5))), 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
k2 = line(basin_lat, basin_Z(:, find(basin_lon == -55.5), find(dens_grid == 27)), 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--');
contour(basin_lat, z_1000, squeeze(basin_gamma_n_z(:, find(basin_lon == -55.5), :))', 'ShowText', 'on');
shading flat; grid off;
xlabel('Latitude °N');
ylabel('Depth');
ylim([-800 0])
xlim([10 35]);
cmocean('balance');
c1 = colorbar;
c1.Label.String = 'mol N m^-^4';
caxis([-2e-8 2e-8]);
title('f) Lat ∇_i_s_o N at 56°W');
lgd = legend([k1 k2],'Winter mixed layer','Thermocline','location','southeast')

%% Plot 2D sections & maps for rates of mixing.


%Remove Kdia and w* values MLD to surface
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

figure
t = tiledlayout(3,2,'TileSpacing','Compact');

% Plot 1: Kiso (m^2 s^-^1) MAP
nexttile(1)
m_proj('mercator', 'lat', [min_lat max_lat],'long', [min_lon max_lon]);
m_pcolor(basin_lon, basin_lat, squeeze(basin_K_iso_z(:,:,find(z_1000 == -200))))
%m_line([-55.5, -55.5], [10, 35],'Color', 'k', 'LineWidth', 2)
m_line(NA_lon_chl, NA_lat_chl, 'color','k','linewi',0.3)
m_hatch(NA_lon_chl, NA_lat_chl,'single',30,5,'color','k'); 
m_grid('fontsize',10);
m_coast('color','k')
set(gcf,'color','w')
set(gca,'fontsize',10)
caxis([500 2500]);
title('a) K_i_s_o at 200m')

% Plot 2: Kiso (m^2 s^-^1) SECTION
nexttile(2)
pcolor(basin_lat, z_1000, squeeze(basin_K_iso_z(:, find(basin_lon == -55.5), :))')
hold on
%line(basin_lat, -squeeze(basin_w_mld_z(:, find(basin_lon == -55.5))), 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
%line(basin_lat, basin_Z(:, find(basin_lon == -55.5), find(dens_grid == 27)), 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--');
contour(basin_lat, z_1000, squeeze(basin_gamma_n_z(:, find(basin_lon == -55.5), :))', 'Color','k','ShowText', 'on');
shading flat; grid off;
xlabel('Latitude °N');
ylabel('Depth');
ylim([-800 0])
xlim([10 35]);
caxis([500 2500]);
cmocean('amp');
c2 = colorbar;
c2.Label.String = 'm^2 s^-^1';
title('b) K_i_s_o at 56°W');
hold on
k1 = line(basin_lat, -squeeze(basin_w_mld_z(:, find(basin_lon == -55.5))), 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
k2 = line(basin_lat, basin_Z(:, find(basin_lon == -55.5), find(dens_grid == 27)), 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--');
lgd = legend([k1 k2],'Winter mixed layer','Thermocline','location','southeast')

% Plot 3: Kdia (m^2 s^-^1) MAP
nexttile(3)
m_proj('mercator', 'lat', [min_lat max_lat],'long', [min_lon max_lon]);
m_pcolor(basin_lon, basin_lat, squeeze(basin_K_dia_z(:,:,find(z_1000 == -200))))
%m_line([-55.5, -55.5], [10, 35],'Color', 'k', 'LineWidth', 2)
m_line(NA_lon_chl, NA_lat_chl, 'color','k','linewi',0.3)
m_hatch(NA_lon_chl, NA_lat_chl,'single',30,5,'color','k'); 
m_grid('fontsize',10);
m_coast('color','k')
set(gcf,'color','w')
set(gca,'fontsize',10)
cmocean('amp');
caxis([0 0.3e-4]);
title('c) K_d_i_a at 200m')

% Plot 4: Kdia (m^2 s^-^1)
nexttile(4)
pcolor(basin_lat, z_1000, squeeze(basin_K_dia_z(:, find(basin_lon == -55.5), :))')
hold on
k1 = line(basin_lat, -squeeze(basin_w_mld_z(:, find(basin_lon == -55.5))), 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
k2 = line(basin_lat, basin_Z(:, find(basin_lon == -55.5), find(dens_grid == 27)), 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--');
contour(basin_lat, z_1000, squeeze(basin_gamma_n_z(:, find(basin_lon == -55.5), :))', 'ShowText', 'on');
shading flat; grid off;
xlabel('Latitude °N');
ylabel('Depth');
ylim([-800 0])
xlim([10 35]);
cmocean('amp');
caxis([0 0.3e-4]);
c3 = colorbar;
c3.Label.String = 'm^2 s^-^1';
title('d) K_d_i_a at 56°W');

% Plot 5: w* (m s^-^1) MAP
nexttile(5)
m_proj('mercator', 'lat', [min_lat max_lat],'long', [min_lon max_lon]);
m_pcolor(basin_lon, basin_lat, squeeze(basin_w_z(:,:,find(z_1000 == -200))))
%m_line([-55.5, -55.5], [10, 35],'Color', 'k', 'LineWidth', 2)
m_line(NA_lon_chl, NA_lat_chl, 'color','k','linewi',0.3)
m_hatch(NA_lon_chl, NA_lat_chl,'single',30,5,'color','k'); 
m_grid('fontsize',10);
m_coast('color','k')
set(gcf,'color','w')
set(gca,'fontsize',10)
cmocean('amp');
caxis([-2e-6 2e-6]);
title('e) w* at 200m')

% Plot 6: w* (m s^-^1) SECTION
nexttile(6)
pcolor(basin_lat, z_1000, squeeze(basin_w_z(:, find(basin_lon == -55.5), :))')
hold on
k1 = line(basin_lat, -squeeze(basin_w_mld_z(:, find(basin_lon == -55.5))), 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
k2 = line(basin_lat, basin_Z(:, find(basin_lon == -55.5), find(dens_grid == 27)), 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--');
contour(basin_lat, z_1000, squeeze(basin_gamma_n_z(:, find(basin_lon == -55.5), :))', 'ShowText', 'on');
shading flat; grid off;
xlabel('Latitude °N');
ylabel('Depth');
ylim([-800 0])
xlim([10 35]);
cmocean('amp');
caxis([-2e-6 2e-6]);
c4 = colorbar;
c4.Label.String = 'm s^-^1';
title('f) w* at 56°W');

%% 2D sections and maps for fluxes

% Remove mixed layer Fdia fluxes
for i = 1:size(basin_F_dia_diff_z, 1)  % Loop over the x-dimension
    for j = 1:size(basin_F_dia_diff_z, 2)  % Loop over the y-dimension
        
        % Find the MLD at the current (x, y) position
        mld = -basin_w_mld_z(i, j);
        
        % Loop over the depth levels (z-dimension)
        for k = 1:size(basin_F_dia_diff_z, 3)
            % If the depth is shallower than the MLD, set the value to NaN
            if z_1000(k) > mld
                basin_F_dia_diff_z(i, j, k) = NaN;
                basin_F_dia_adv_z(i,j,k) = NaN;
            end
        end
        
    end
end

figure
t = tiledlayout(4,2,'TileSpacing','Compact');

% Plot 1: Longitudinal Fiso (mol m^2 yr^-^1) MAP
nexttile(1)
m_proj('mercator', 'lat', [min_lat max_lat],'long', [min_lon max_lon]);
m_pcolor(basin_lon, basin_lat, squeeze(basin_F_iso_lon_z(:,:,find(z_1000 == -200))))
%m_line([-55.5, -55.5], [10, 35],'Color', 'k', 'LineWidth', 2)
m_line(NA_lon_chl, NA_lat_chl, 'color','k','linewi',0.5)
m_hatch(NA_lon_chl, NA_lat_chl,'single',30,5,'color','k'); 
m_grid('fontsize',10);
m_coast('color','k')
set(gcf,'color','w')
set(gca,'fontsize',10)
cmocean('balance')
caxis([-250 250])
title('a) Longitudinal F_i_s_o at 200m')

% Plot 2: Longitudinal Fiso (mol m^2 yr^-^1) SECTION
nexttile(2)
pcolor(basin_lat, z_1000, squeeze(basin_F_iso_lon_z(:, find(basin_lon == -55.5), :))')
hold on
k1 = line(basin_lat, -squeeze(basin_w_mld_z(:, find(basin_lon == -55.5))), 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
k2 = line(basin_lat, basin_Z(:, find(basin_lon == -55.5), find(dens_grid == 27)), 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--');
contour(basin_lat, z_1000, squeeze(basin_gamma_n_z(:, find(basin_lon == -55.5), :))', 'Color','k','ShowText', 'on');
shading flat; grid off;
xlabel('Latitude °N');
ylabel('Depth');
ylim([-800 0])
xlim([10 35]);
cmocean('balance')
caxis([-250 250])
c = colorbar; 
c.Label.String = 'mol m^2 yr^-^1';
title('b) Longitudinal F_i_s_o at 56°W')
lgd = legend([k1 k2],'Winter mixed layer','Thermocline','location','southeast')

% Plot 3: Latitudinal Fiso (mol m^2 yr^-^1) MAP
nexttile(3)
m_proj('mercator', 'lat', [min_lat max_lat],'long', [min_lon max_lon]);
m_pcolor(basin_lon, basin_lat, squeeze(basin_F_iso_lat_z(:,:,find(z_1000 == -200))))
%m_line([-55.5, -55.5], [10, 35],'Color', 'k', 'LineWidth', 2)
m_line(NA_lon_chl, NA_lat_chl, 'color','k','linewi',0.5)
m_hatch(NA_lon_chl, NA_lat_chl,'single',30,5,'color','k'); 
m_grid('fontsize',10);
m_coast('color','k')
set(gcf,'color','w')
set(gca,'fontsize',10)
cmocean('balance')
caxis([-250 250])
title('c) Latitudinal F_i_s_o at 200m')

% Plot 4: Latitudinal Fiso (mol m^2 yr^-^1) SECTION
nexttile(4)
pcolor(basin_lat, z_1000, squeeze(basin_F_iso_lat_z(:, find(basin_lon == -55.5), :))')
hold on
k1 = line(basin_lat, -squeeze(basin_w_mld_z(:, find(basin_lon == -55.5))), 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
k2 = line(basin_lat, basin_Z(:, find(basin_lon == -55.5), find(dens_grid == 27)), 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--');
contour(basin_lat, z_1000, squeeze(basin_gamma_n_z(:, find(basin_lon == -55.5), :))', 'Color','k', 'ShowText', 'on');
shading flat; grid off;
xlabel('Latitude °N');
ylabel('Depth');
ylim([-800 0])
xlim([10 35]);
cmocean('balance');
caxis([-250 250]);
c = colorbar;
c.Label.String = 'mol m^2 yr^-^1'
title('d) Latitudinal F_i_s_o at 56°W');

% Plot 5:Fdia (mol m^2 yr^-^1) MAP
nexttile(5)
m_proj('mercator', 'lat', [min_lat max_lat],'long', [min_lon max_lon]);
m_pcolor(basin_lon, basin_lat, squeeze(basin_F_dia_diff_z(:,:,find(z_1000 == -200))))
%m_line([-55.5, -55.5], [10, 35],'Color', 'k', 'LineWidth', 2)
m_line(NA_lon_chl, NA_lat_chl, 'color','k','linewi',0.5)
m_hatch(NA_lon_chl, NA_lat_chl,'single',30,5,'color','k'); 
m_grid('fontsize',10);
m_coast('color','k')
set(gcf,'color','w')
set(gca,'fontsize',10)
cmocean('amp');
cmocean('balance');
caxis([-0.07 0.07]);
title('e) F_d_i_a at 200m')

% Plot 6: Fdia (mol m^2 yr^-^1) SECTION
nexttile(6)
pcolor(basin_lat, z_1000, squeeze(basin_F_dia_diff_z(:, find(basin_lon == -55.5), :))')
hold on
k1 = line(basin_lat, -squeeze(basin_w_mld_z(:, find(basin_lon == -55.5))), 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
k2 = line(basin_lat, basin_Z(:, find(basin_lon == -55.5), find(dens_grid == 27)), 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--');
contour(basin_lat, z_1000, squeeze(basin_gamma_n_z(:, find(basin_lon == -55.5), :))', 'Color','k', 'ShowText', 'on');
shading flat; grid off;
xlabel('Latitude °N');
ylabel('Depth');
ylim([-800 0])
xlim([10 35]);
cmocean('balance');
caxis([-0.07 0.07]);
c=colorbar;
c.Label.String = 'mol m^2 yr^-^1'
title('f) F_d_i_a at 56°W');

% Plot 7: w*N (mol m^2 yr^-^1) MAP
nexttile(7)
m_proj('mercator', 'lat', [min_lat max_lat],'long', [min_lon max_lon]);
m_pcolor(basin_lon, basin_lat, squeeze(basin_F_dia_adv_z(:,:,find(z_1000 == -200))))
%m_line([-55.5, -55.5], [10, 35],'Color', 'k', 'LineWidth', 2)
m_line(NA_lon_chl, NA_lat_chl, 'color','k','linewi',0.5)
m_hatch(NA_lon_chl, NA_lat_chl,'single',30,5,'color','k'); 
m_grid('fontsize',10);
m_coast('color','k')
set(gcf,'color','w')
set(gca,'fontsize',10)
cmocean('balance')
caxis([-0.2 0.2])
title('g) w*N at 200m')

% Plot 8: w*N (mol m^2 yr^-^1) SECTION
nexttile(8)
pcolor(basin_lat, z_1000, squeeze(basin_F_dia_adv_z(:, find(basin_lon == -55.5), :))')
hold on
k1 = line(basin_lat, -squeeze(basin_w_mld_z(:, find(basin_lon == -55.5))), 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
k2 = line(basin_lat, basin_Z(:, find(basin_lon == -55.5), find(dens_grid == 27)), 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--');
contour(basin_lat, z_1000, squeeze(basin_gamma_n_z(:, find(basin_lon == -55.5), :))', 'Color','k','ShowText', 'on');
shading flat; grid off;
xlabel('Latitude °N');
ylabel('Depth');
ylim([-800 0])
xlim([10 35]);
cmocean('balance');
caxis([-0.2 0.2]);
c = colorbar
c.Label.String = 'mol m^2 yr^-^1'
title('h) w*N at 56°W');

%% 2D maps of supply (mol m-2 yr-1) to wml and upper thermocline

%eddy stirring
figure
t = tiledlayout(2,2,'TileSpacing','Compact');
nexttile(1)
m_proj('mercator', 'lat', [min_lat max_lat],'long', [min_lon max_lon]);
m_pcolor(basin_lon, basin_lat, basin_wml_F_iso)%_area)
hold on
m_line(NA_lon_chl, NA_lat_chl, 'color','k','linewi',0.3)
m_hatch(NA_lon_chl, NA_lat_chl,'single',30,5,'color','k'); 
m_grid('fontsize',6);
m_coast('color','k')
set(gcf,'color','w')
cmocean('balance')
caxis([-0.17 0.17])
c = colorbar;
c.Label.String = 'mol m^-^2 yr^-^1';
title('a) Eddy stirring supply','Units', 'normalized', 'Position', [0.5, 1.1, 1.2])
set(gca,'Fontsize',7)
%microscale turbulence
nexttile(2)
m_proj('mercator', 'lat', [min_lat max_lat],'long', [min_lon max_lon]);
m_pcolor(basin_lon, basin_lat, (basin_wml_F_dia_diff + basin_wml_F_dia_adv))%_area + basin_wml_F_dia_adv_area))
hold on
m_line(NA_lon_chl, NA_lat_chl, 'color','k','linewi',0.3)
m_hatch(NA_lon_chl, NA_lat_chl,'single',30,5,'color','k'); 
m_grid('fontsize',6);
m_coast('color','k');
set(gcf,'color','w')
cmocean('balance')
caxis([-0.17 0.17])
c = colorbar;
c.Label.String = 'mol m^-^2 yr^-^1';
title('b) Turbulent mixing supply','Units', 'normalized', 'Position', [0.5, 1.1, 1.2])
set(gca,'Fontsize',7)
% subplot(2,2,3)
nexttile(3)
m_proj('mercator', 'lat', [min_lat max_lat],'long', [min_lon max_lon]);
m_pcolor(basin_lon, basin_lat, basin_therm_F_iso)%_area)
hold on
m_line(NA_lon_chl, NA_lat_chl, 'color','k','linewi',0.3)
m_hatch(NA_lon_chl, NA_lat_chl,'single',30,5,'color','k'); 
m_grid('fontslatze',6);
m_coast('color','k');
set(gcf,'color','w')
cmocean('balance')
caxis([-0.6 0.6])%caxis([-40 40])
c = colorbar;
c.Label.String = 'mol m^-^2 yr^-^1';
set(gca,'Fontsize',7)
% Define the tick positions and labels for the colorbar
tickPositions = -0.6:0.2:0.6;
c.Ticks = tickPositions;
c.TickLabels = string(tickPositions);
% subplot(2,2,4)
nexttile(4)
m_proj('mercator', 'lat', [min_lat max_lat],'long', [min_lon max_lon]);
m_pcolor(basin_lon, basin_lat, (basin_therm_F_dia_diff + basin_therm_F_dia_adv))%_area + basin_therm_F_dia_adv_area))
hold on
m_line(NA_lon_chl, NA_lat_chl, 'color','k','linewi',0.3)
m_hatch(NA_lon_chl, NA_lat_chl,'single',30,5,'color','k'); 
m_grid('fontsize',6);
m_coast('color','k');
set(gcf,'color','w')
cmocean('balance')
caxis([-0.17 0.17])%caxis([-40 40])
c = colorbar;
c.Label.String = 'mol m^-^2 yr^-^1';
set(gca,'Fontsize',7)

%explore values of Fiso thermocline supply maps
% section = (basin_wml_F_iso.*basin_NA_mask)
% 
% lon = [min_lon:max_lon]
% lat = [min_lat:max_lat]
% 
% mAx = max(section(:));
% mIn= min(section(:));
% abs_section = abs(section);
% mIn_abs = min(abs_section(:));
% 
% [mAx_lat, mAx_lon] = find(ismember(section, max(section(:))));
% [mIn_lat, mIn_lon] = find(ismember(section, min(section(:))));
% [mIn_abs_lat, mIn_abs_lon] = find(ismember(abs_section, min(abs_section(:))));
% 
% stats(1,1) = mAx
% stats(1,2) = lat(mAx_lat)
% stats(1,3) = lon(mAx_lon)
% stats(2,1) = mIn
% stats(2,2) = lat(mIn_lat)
% stats(2,3) = lon(mIn_lon)
% stats(3,1) = mIn_abs
% 
% %-0.0746   14.5000  -44.5000 WML loss
% % -0.0701   14.5000  -44.5000
% 
% 
% %% Basin 3D slices (N, mixing, fluxes)
% 
% %cut depth of basin area from 1000m to 400m:
% basin_N_3D = basin_N_z(:,:,[1:find(z_1000 == -400)]);
% basin_N_grad_dia_3D = basin_N_grad_dia_z(:,:,[1:find(z_1000 == -400)]);
% basin_N_grad_iso_lon_3D = basin_N_grad_iso_lon_z(:,:,[1:find(z_1000 == -400)]);
% basin_N_grad_iso_lat_3D = basin_N_grad_iso_lat_z(:,:,[1:find(z_1000 == -400)]);
% basin_K_dia_3D = basin_K_dia_z(:,:,[1:find(z_1000 == -400)]);
% basin_K_iso_3D = basin_K_iso_z(:,:,[1:find(z_1000 == -400)]);
% basin_w_3D = basin_w_z(:,:,[1:find(z_1000 == -400)]);
% basin_F_dia_diff_3D = basin_F_dia_diff_z(:,:,[1:find(z_1000 == -400)]);
% basin_F_dia_adv_3D = basin_F_dia_adv_z(:,:,[1:find(z_1000 == -400)]);
% basin_F_iso_lon_3D = basin_F_iso_lon_z(:,:,[1:find(z_1000 == -400)]);
% basin_F_iso_lat_3D = basin_F_iso_lat_z(:,:,[1:find(z_1000 == -400)]);
% basin_conv_F_iso_3D = basin_conv_F_iso_z(:,:,[1:find(z_1000 == -400)]);
% basin_conv_F_dia_diff_3D = basin_conv_F_dia_diff_z(:,:,[1:find(z_1000 == -400)]);
% basin_conv_F_dia_adv_3D = basin_conv_F_dia_adv_z(:,:,[1:find(z_1000 == -400)]);
% basin_gamma_n_z_3D = basin_gamma_n_z(:,:,[1:find(z_1000 == -400)]);
% 
% %set up slices
% [x,y,z] = meshgrid(basin_lon,basin_lat,double(z_1000([1:find(z_1000 == -400)])));
% 
% basin = input('Enter basin: ')
% 
% switch basin
% 
%     case 'NA'
%         xslice = [-19.5 -40.5 -60.5];                
%         yslice = 0;%40;     
%         zslice = [-400 -200]; 
%         x_lower_lim = -80;
%         x_upper_lim = -19.5;
%         y_lower_lim = 8;
%         y_upper_lim = 40;
%         gyre_lon = NA_lon_chl;
%         gyre_lat = NA_lat_chl;
%     case 'SA'
%         xslice = [4.5 -10.5 -30.5];
%         yslice = -0.5;
%         zslice = [-400 -200];
%         x_lower_lim = -40
%         x_upper_lim = 5.5
%         y_lower_lim = -40
%         y_upper_lim = -0.5
%         gyre_lon = SA_lon_chl;
%         gyre_lat = SA_lat_chl;
%     case 'I'
%         xslice = [110 85 60];
%         yslice = -9.5;
%         zslice = [-400 -200];
%         x_lower_lim = 45
%         x_upper_lim = 110
%         y_lower_lim = -40
%         y_upper_lim = -9.5
%         gyre_lon = I_lon_chl;
%         gyre_lat = I_lat_chl;
%     case 'NP'
%         xslice = [230 200 160];
%         yslice = 35;
%         zslice = [-400 -200];
%         x_lower_lim = 115
%         x_upper_lim = 230
%         y_lower_lim = 0.5
%         y_upper_lim = 35
%         gyre_lon = [NP1_lon_chl; (NP2_lon_chl + 360)]
%         gyre_lat = [NP1_lat_chl; NP2_lat_chl]
%      case 'SP'
%         xslice = [280 230 190];
%         yslice = -0.5;
%         zslice = [-400 -200];
%         x_lower_lim = 140;
%         x_upper_lim = 280;
%         y_lower_lim = -40.5;
%         y_upper_lim = -0.5;
%         SP1_lat_chl_fix = flip([SP1_lat_chl([32:end]);SP1_lat_chl([1:31])]);
%         SP1_lon_chl_fix = flip([SP1_lon_chl([32:end]);SP1_lon_chl([1:31])]);
%         gyre_lon = [SP1_lon_chl_fix; (SP2_lon_chl + 360)];
%         gyre_lat = [SP1_lat_chl_fix; SP2_lat_chl];
% end
% 
% %% Plot basin N distribution 
% 
% figure
% subplot(2,2,1)
% s = slice(x,y,z,basin_N_3D,xslice,yslice,zslice)
% set(s,'LineWidth',0.0001)
% xlabel('Longitude °E')
% ylabel('Latitude °N')
% zlabel('Depth (m)')
% xlim([x_lower_lim x_upper_lim])
% ylim([y_lower_lim y_upper_lim])
% shading flat
% hold on
% fill3([gyre_lon],[gyre_lat],[repmat(-200,length(gyre_lon),1)],'k','FaceColor','none','EdgeColor','k','LineWidth',1)
% fill3([gyre_lon],[gyre_lat],[repmat(-400,length(gyre_lon),1)],'k','FaceColor','none','EdgeColor','k','LineWidth',1)
% contourslice(x,y,z,basin_gamma_n_z_3D,xslice,[],[])%,'ShowText','on')
% hold off
% cmocean('algae')
% a = colorbar
% a.Label.String = 'mol N m^-^3'
% title('a) N')
% view([-15.05172413793108,15.585992217898843])
% caxis([0 0.03])
% 
% subplot(2,2,2)
% s = slice(x,y,z,basin_N_grad_iso_lon_3D,xslice,yslice,zslice)
% set(s,'LineWidth',0.0001)
% xlabel('Longitude °E')
% ylabel('Latitude °N')
% zlabel('Depth (m)')
% xlim([x_lower_lim x_upper_lim])
% ylim([y_lower_lim y_upper_lim])
% shading flat
% hold on
% fill3([gyre_lon],[gyre_lat],[repmat(-200,length(gyre_lon),1)],'k','FaceColor','none','EdgeColor','k','LineWidth',1)
% fill3([gyre_lon],[gyre_lat],[repmat(-400,length(gyre_lon),1)],'k','FaceColor','none','EdgeColor','k','LineWidth',1)
% contourslice(x,y,z,basin_gamma_n_z_3D,xslice,[],[])%,'ShowText','on')
% hold off
% grid on
% cmocean('balance')
% caxis([-2e-8 2e-8])
% c = colorbar
% cmocean('balance')
% c.Label.String = 'mol N m^-^4'
% title('b) Longitudinal  ∇_i_s_o N')
% view([-15.05172413793108,15.585992217898843])
% 
% subplot(2,2,3)
% s = slice(x,y,z,basin_N_grad_iso_lat_3D,xslice,yslice,zslice)
% set(s,'LineWidth',0.0001)
% xlabel('Longitude °E')
% ylabel('Latitude °N')
% zlabel('Depth (m)')
% xlim([x_lower_lim x_upper_lim])
% ylim([y_lower_lim y_upper_lim])
% shading flat
% hold on
% fill3([gyre_lon],[gyre_lat],[repmat(-200,length(gyre_lon),1)],'k','FaceColor','none','EdgeColor','k','LineWidth',1)
% fill3([gyre_lon],[gyre_lat],[repmat(-400,length(gyre_lon),1)],'k','FaceColor','none','EdgeColor','k','LineWidth',1)
% contourslice(x,y,z,basin_gamma_n_z_3D,xslice,[],[])
% hold off
% grid on
% cmocean('balance')
% caxis([-2e-8 2e-8])
% c = colorbar
% c.Label.String = 'mol N m^-^4'
% title('c) Latitudinal ∇_i_s_o N')
% view([-15.05172413793108,15.585992217898843])
% 
% subplot(2,2,4)
% s = slice(x,y,z,basin_N_grad_dia_3D,xslice,yslice,zslice)
% set(s,'LineWidth',0.0001)
% xlabel('Longitude °E')
% ylabel('Latitude °N')
% zlabel('Depth (m)')
% xlim([x_lower_lim x_upper_lim])
% ylim([y_lower_lim y_upper_lim])
% shading flat
% hold on
% fill3([gyre_lon],[gyre_lat],[repmat(-200,length(gyre_lon),1)],'k','FaceColor','none','EdgeColor','k','LineWidth',1)
% fill3([gyre_lon],[gyre_lat],[repmat(-400,length(gyre_lon),1)],'k','FaceColor','none','EdgeColor','k','LineWidth',1)
% contourslice(x,y,z,basin_gamma_n_z_3D,xslice,[],[])
% cmocean('balance')
% c = colorbar
% caxis([-1e-4 1e-4])
% c.Label.String = 'mol N m^-^4'
% title('d) dN/dz*')
% view([-15.05172413793108,15.585992217898843])
% 
% %% Plot rates of mixing
% 
% %First, remove Kdia and w* values MLD to surface
% % Assuming:
% % basin_K_dia_3D is the 3D matrix of size (37x71x33) corresponding to (x, y, z)
% % basin_w_mld_z is the 2D matrix of size (37x71) corresponding to (x, y)
% % z is a 1D array of size 33 corresponding to the depth levels
% z_1D = squeeze(z(1,1,:))
% 
% % Loop over the x and y dimensions
% for i = 1:size(basin_K_dia_3D, 1)  % Loop over the x-dimension
%     for j = 1:size(basin_K_dia_3D, 2)  % Loop over the y-dimension
% 
%         % Find the MLD at the current (x, y) position
%         mld = -basin_w_mld_z(i, j);
% 
%         % Loop over the depth levels (z-dimension)
%         for k = 1:size(basin_K_dia_3D, 3)
%             % If the depth is shallower than the MLD, set the value to NaN
%             if z_1D(k) > mld
%                 basin_K_dia_3D(i, j, k) = NaN;
%                 basin_w_3D(i,j,k) = NaN;
%             end
%         end
% 
%     end
% end
% 
% 
% %Rates of Mixing
% figure
% t = tiledlayout(3,1,'TileSpacing','Compact');
% nexttile(1)
% s = slice(x,y,z,basin_K_iso_3D,xslice,yslice,zslice)
% set(s,'LineWidth',0.0001)
% xlabel('Longitude °E')
% ylabel('Latitude °N')
% zlabel('Depth (m)')
% xlim([x_lower_lim x_upper_lim])
% ylim([y_lower_lim y_upper_lim])
% shading flat
% hold on
% fill3([gyre_lon],[gyre_lat],[repmat(-200,length(gyre_lon),1)],'k','FaceColor','none','EdgeColor','k','LineWidth',1)
% fill3([gyre_lon],[gyre_lat],[repmat(-400,length(gyre_lon),1)],'k','FaceColor','none','EdgeColor','k','LineWidth',1)
% s = contourslice(x,y,z,basin_gamma_n_z_3D,xslice,[],[])%,'ShowText','on')
% set(s,'EdgeColor','k')
% hold off
% cmocean('amp')
% caxis([0 2000])
% c = colorbar
% c.Label.String = 'm^2 s^-^1'
% title('a) K_i_s_o')
% view([-15.05172413793108,15.585992217898843])
% 
% nexttile(2)
% s = slice(x,y,z,basin_K_dia_3D,xslice,yslice,zslice)
% set(s,'LineWidth',0.0001)
% xlabel('Longitude °E')
% ylabel('Latitude °N')
% zlabel('Depth (m)')
% xlim([x_lower_lim x_upper_lim])
% ylim([y_lower_lim y_upper_lim])
% shading flat
% hold on
% fill3([gyre_lon],[gyre_lat],[repmat(-200,length(gyre_lon),1)],'k','FaceColor','none','EdgeColor','k','LineWidth',1)
% fill3([gyre_lon],[gyre_lat],[repmat(-400,length(gyre_lon),1)],'k','FaceColor','none','EdgeColor','k','LineWidth',1)
% s = contourslice(x,y,z,basin_gamma_n_z_3D,xslice,[],[])%,'ShowText','on')
% set(s,'EdgeColor','k')
% hold off
% cmocean('amp')
% c = colorbar
% caxis([0 0.3e-4])
% c.Label.String = 'm^2 s^-^1'
% title('b) K_d_i_a')
% view([-15.05172413793108,15.585992217898843])
% 
% nexttile(3)
% s = slice(x,y,z,basin_w_3D,xslice,yslice,zslice)
% set(s,'LineWidth',0.0001)
% xlabel('Longitude °E')
% ylabel('Latitude °N')
% zlabel('Depth (m)')
% xlim([x_lower_lim x_upper_lim])
% ylim([y_lower_lim y_upper_lim])
% shading flat
% hold on
% fill3([gyre_lon],[gyre_lat],[repmat(-200,length(gyre_lon),1)],'k','FaceColor','none','EdgeColor','k','LineWidth',1)
% fill3([gyre_lon],[gyre_lat],[repmat(-400,length(gyre_lon),1)],'k','FaceColor','none','EdgeColor','k','LineWidth',1)
% s = contourslice(x,y,z,basin_gamma_n_z_3D,xslice,[],[])%,'ShowText','on')
% set(s,'EdgeColor','k')
% hold off
% cmocean('balance')
% c = colorbar
% caxis([-2e-6 2e-6])
% c.Label.String = 'm s^-^1'
% title('c) w*')
% view([-15.05172413793108,15.585992217898843])
% 
% %Calculate max and min of w* in 3D section
% 
% my_stats_slice(basin_w_3D, basin_gyre_mask)
% %max of 7.1539e-06 at 32.5N 38.5W and 20m
% %max of -8.7293e-06 at 22.5N 36.5W and 5m
% 
% %% Plot basin fluxes
% 
% 
% %First, remove F_dia_diff and F_dia_adv above the MLD
% z_1D = squeeze(z(1,1,:))
% % Loop over the x and y dimensions
% for i = 1:size(basin_F_dia_diff_3D, 1)  % Loop over the x-dimension
%     for j = 1:size(basin_F_dia_diff_3D, 2)  % Loop over the y-dimension
% 
%         % Find the MLD at the current (x, y) position
%         mld = -basin_w_mld_z(i, j);
% 
%         % Loop over the depth levels (z-dimension)
%         for k = 1:size(basin_F_dia_diff_3D, 3)
%             % If the depth is shallower than the MLD, set the value to NaN
%             if z_1D(k) > mld
%                 basin_F_dia_diff_3D(i, j, k) = NaN;
%                 basin_F_dia_adv_3D(i,j,k) = NaN;
%             end
%         end
% 
%     end
% end
% 
% figure
% t = tiledlayout(2,2,'TileSpacing','Compact')
% nexttile
% s = slice(x,y,z,basin_F_iso_lon_3D,xslice,yslice,zslice)
% set(s,'LineWidth',0.0001)
% xlabel('Longitude °E')
% ylabel('Latitude °N')
% zlabel('Depth (m)')
% xlim([x_lower_lim x_upper_lim])
% ylim([y_lower_lim y_upper_lim])
% shading flat
% hold on
% fill3([gyre_lon],[gyre_lat],[repmat(-200,length(gyre_lon),1)],'k','FaceColor','none','EdgeColor','k','LineWidth',1)
% fill3([gyre_lon],[gyre_lat],[repmat(-400,length(gyre_lon),1)],'k','FaceColor','none','EdgeColor','k','LineWidth',1)
% s = contourslice(x,y,z,basin_gamma_n_z_3D,xslice,[],[])%,'ShowText','on')
% set(s,'EdgeColor','k')
% hold off
% cmocean('amp')
% c = colorbar
% cmocean('balance')
% caxis([-250 250])
% c.Label.String = 'mol m^2 yr^-^1'
% title('a) Longitudinal F_i_s_o')
% view([-15.05172413793108,15.585992217898843])
% 
% nexttile
% s = slice(x,y,z,basin_F_iso_lat_3D,xslice,yslice,zslice)
% set(s,'LineWidth',0.0001)
% xlabel('Longitude °E')
% ylabel('Latitude °N')
% zlabel('Depth (m)')
% xlim([x_lower_lim x_upper_lim])
% ylim([y_lower_lim y_upper_lim])
% shading flat
% hold on
% fill3([gyre_lon],[gyre_lat],[repmat(-200,length(gyre_lon),1)],'k','FaceColor','none','EdgeColor','k','LineWidth',1)
% fill3([gyre_lon],[gyre_lat],[repmat(-400,length(gyre_lon),1)],'k','FaceColor','none','EdgeColor','k','LineWidth',1)
% s = contourslice(x,y,z,basin_gamma_n_z_3D,xslice,[],[])%,'ShowText','on')
% set(s,'EdgeColor','k')
% cmocean('amp')
% c = colorbar
% cmocean('balance')
% caxis([-250 250])
% c.Label.String = 'mol m^2 yr^-^1'
% title('b) Latitudinal F_i_s_o')
% view([-15.05172413793108,15.585992217898843])
% 
% %Nitrate gradients
% nexttile
% s = slice(x,y,z,basin_F_dia_diff_3D,xslice,yslice,zslice)
% set(s,'LineWidth',0.0001)
% xlabel('Longitude °E')
% ylabel('Latitude °N')
% zlabel('Depth (m)')
% xlim([x_lower_lim x_upper_lim])
% ylim([y_lower_lim y_upper_lim])
% shading flat
% hold on
% fill3([gyre_lon],[gyre_lat],[repmat(-200,length(gyre_lon),1)],'k','FaceColor','none','EdgeColor','k','LineWidth',1)
% fill3([gyre_lon],[gyre_lat],[repmat(-400,length(gyre_lon),1)],'k','FaceColor','none','EdgeColor','k','LineWidth',1)
% contourslice(x,y,z,basin_gamma_n_z_3D,xslice,[],[])%,'ShowText','on')
% hold off
% grid on
% cmocean('balance')
% caxis([-0.07 0.07])
% c = colorbar
% cmocean('balance')
% c.Label.String = 'mol m^2 yr^-^1'
% title('c) F_d_i_a')
% view([-15.05172413793108,15.585992217898843])
% 
% nexttile
% s = slice(x,y,z,basin_F_dia_adv_3D,xslice,yslice,zslice)
% set(s,'LineWidth',0.0001)
% xlabel('Longitude °E')
% ylabel('Latitude °N')
% zlabel('Depth (m)')
% xlim([x_lower_lim x_upper_lim])
% ylim([y_lower_lim y_upper_lim])
% shading flat
% hold on
% fill3([gyre_lon],[gyre_lat],[repmat(-200,length(gyre_lon),1)],'k','FaceColor','none','EdgeColor','k','LineWidth',1)
% fill3([gyre_lon],[gyre_lat],[repmat(-400,length(gyre_lon),1)],'k','FaceColor','none','EdgeColor','k','LineWidth',1)
% contourslice(x,y,z,basin_gamma_n_z_3D,xslice,[],[])%,'ShowText','on')
% hold off
% grid on
% cmocean('balance')
% caxis([-0.2 0.2])
% c = colorbar
% c.Label.String = 'mol m^2 yr^-^1'
% title('d) w*N')
% view([-15.05172413793108,15.585992217898843])
% 
% my_stats_slice(basin_F_dia_adv_3D, basin_gyre_mask)


