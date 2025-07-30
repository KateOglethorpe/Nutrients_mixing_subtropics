% Kate Oglethorpe
% May 2025

% This script examines the biological definition for the subtropical gyres
% used in the paper: chl-a < 0.1 mg -3. 

close all; clear


%% Calculate mean 2021 chl-a (mg m-3) concentration 

% Load monthly 2021 chl-a from:
% https://neo.gsfc.nasa.gov/view.php?datasetId=MY1DMM_CHLORA&year=2021
% Couldn't find annual data

chl_01 = load('Jan_2021.txt')   
chl_02 = load('Feb_2021.txt');
chl_03 = load('Mar_2021.txt')
chl_04 = load('Apr_2021.txt')
chl_05 = load('May_2021.txt')
chl_06 = load('Jun_2021.txt')
chl_07 = load('Jul_2021.txt')
chl_08 = load('Aug_2021.txt')
chl_09 = load('Sep_2021.txt')
chl_10 = load('Oct_2021.txt')
chl_11 = load('Nov_2021.txt')
chl_12 = load('Dec_2021.txt')

% Calculate annual mean & save
chl = nan(180,360);
for lat = 1:180
    for lon = 1:360
        chl(lat,lon) = mean([chl_01(lat,lon), chl_02(lat,lon), chl_03(lat,lon), ...
            chl_04(lat,lon), chl_05(lat,lon), chl_06(lat,lon), chl_07(lat,lon),...
            chl_08(lat,lon), chl_09(lat,lon), chl_10(lat,lon),chl_11(lat,lon),...
            chl_12(lat,lon)], 'omitnan');
    end
end
save('chl.mat','chl')

% Plot 2021 mean and overlay hand-drawn chl-a that follow the chl-a = 0.1
% mg m-3 contour as best as possible. 
chl(chl == 99999) = NaN;
chl = flipud(chl);
figure
m_proj('miller','lat',[-60 60],'lon',[-180 180]);
m_pcolor(glob_lon, glob_lat, chl)
hold on
m_line(NA_lon_chl, NA_lat_chl, 'color','k','linewi',0.3)
m_hatch(NA_lon_chl, NA_lat_chl,'single',30,6,'color','k'); 
m_line(SA_lon_chl, SA_lat_chl, 'color','k','linewi',0.3)
m_hatch(SA_lon_chl, SA_lat_chl,'single',30,6,'color','k'); 
m_line(I_lon_chl, I_lat_chl, 'color','k','linewi',0.3)
m_hatch(I_lon_chl, I_lat_chl,'single',30,6,'color','k'); 
m_line(NP1_lon_chl, NP1_lat_chl, 'color','k','linewi',0.3)
m_hatch(NP1_lon_chl, NP1_lat_chl,'single',30,6,'color','k'); 
m_line(SP1_lon_chl, SP1_lat_chl, 'color','k','linewi',0.3)
m_hatch(SP1_lon_chl, SP1_lat_chl,'single',30,6,'color','k'); 
m_line(NP2_lon_chl, NP2_lat_chl, 'color','k','linewi',0.3)
m_hatch(NP2_lon_chl, NP2_lat_chl,'single',30,6,'color','k'); 
m_line(SP2_lon_chl, SP2_lat_chl, 'color','k','linewi',0.3)
m_hatch(SP2_lon_chl, SP2_lat_chl,'single',30,6,'color','k'); 
m_line([x_lower_lim; x_upper_lim; x_upper_lim; x_lower_lim; x_lower_lim],[y_upper_lim; y_upper_lim; y_lower_lim; y_lower_lim; y_upper_lim],'color','r','linewi',1)
m_grid('fontsize',9);
m_coast('color','k')
set(gcf,'color','w')
c = colorbar
cmocean('algae')
caxis([0 0.3])
c.Label.String = 'Chl-a concentration (mg m^-^3)'
set(gca,'Fontsize',12)

%% Compare chl-a = 0.1 mg m-3 contour for 2021 to other years

% Load 1997-2002 SeaWiFiS data
file = '/Users/ko389/Library/CloudStorage/OneDrive-UniversityofCambridge/MSc_paper/Code/data/chl/SeaWiFS_chl.nc';
lat  = ncread(file, 'lat');
lon  = ncread(file, 'lon');
time = ncread(file, 'time');
data = ncread(file, 'data');

% Apply scale factor and mask missing values
scale_factor = ncreadatt(file, 'data', 'scale_factor');
missing_value = ncreadatt(file, 'data', 'missing_value');
data = double(data);
data(data >30) = NaN;
%data = data * scale_factor;

% Convert time to datetime
time_ref = datetime(1997,1,1,1,0,0);
time = time_ref + hours(time);

% Identify unique years
years = unique(year(time));

% Preallocate annual mean matrix: lon x lat x year
annual_mean = NaN(length(lon), length(lat),length(years));

% Loop through each year
for i = 1:length(years)
    % Find indices for current year
    idx = year(time) == years(i);
    
    % Calculate mean over time dimension (3rd dim)
    annual_mean(:,:,i) = mean(data(:,:,idx), 3, 'omitnan');
end

%Change dimensions from 360x180x6 to 180x360x6
annual_mean = permute(annual_mean, [2 1 3]);

%Flip in latitude direction
annual_mean = flip(annual_mean,1);
%Reorder longitude
annual_mean= cat(2, annual_mean(:,181:360,:), annual_mean(:,1:180,:));

% Plot 0.1 contours for 1997-2002 chl data
figure
m_proj('miller');
m_pcolor(glob_lon, glob_lat, double(chl<0.1))
hold on;
m_contour(glob_lon, glob_lat, annual_mean(:,:,years==1997),[0.1 0.1],'r','LineWidth',0.5)
m_contour(glob_lon, glob_lat, annual_mean(:,:,years==1998),[0.1 0.1],'r','LineWidth',0.5)
m_contour(glob_lon, glob_lat, annual_mean(:,:,years==1999),[0.1 0.1],'r','LineWidth',0.5)
m_contour(glob_lon, glob_lat, annual_mean(:,:,years==2000),[0.1 0.1],'r','LineWidth',0.5)
m_contour(glob_lon, glob_lat, annual_mean(:,:,years==2001),[0.1 0.1],'r','LineWidth',0.5)
m_contour(glob_lon, glob_lat, annual_mean(:,:,years==2002),[0.1 0.1],'r','LineWidth',0.5)
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
%c = colorbar
%cmocean('algae')
%caxis([0 0.1])
%c.Label.String = 'Chl-a concentration (mg m^-^3)'
set(gca,'Fontsize',9)

% Load 2002-2022 chl data from ESA 4km resolution gridded product:
% https://data.ceda.ac.uk/neodc/esacci/ocean_colour/data/v6.0-release/geographic/netcdf/chlor_a/annual/v6.0
% File path
ncfile_2022 = 'ESACCI-OC-L3S-CHLOR_A-MERGED-1Y_YEARLY_4km_GEO_PML_OCx-2022-fv6.0.nc';    %2022
ncfile_2021 = 'ESACCI-OC-L3S-CHLOR_A-MERGED-1Y_YEARLY_4km_GEO_PML_OCx-2021-fv6.0.nc';
ncfile_2020 = 'ESACCI-OC-L3S-CHLOR_A-MERGED-1Y_YEARLY_4km_GEO_PML_OCx-2020-fv6.0.nc';
ncfile_2019 = 'ESACCI-OC-L3S-CHLOR_A-MERGED-1Y_YEARLY_4km_GEO_PML_OCx-2019-fv6.0.nc';
ncfile_2018 = 'ESACCI-OC-L3S-CHLOR_A-MERGED-1Y_YEARLY_4km_GEO_PML_OCx-2018-fv6.0.nc';
ncfile_2017 = 'ESACCI-OC-L3S-CHLOR_A-MERGED-1Y_YEARLY_4km_GEO_PML_OCx-2017-fv6.0.nc';
ncfile_2016 = 'ESACCI-OC-L3S-CHLOR_A-MERGED-1Y_YEARLY_4km_GEO_PML_OCx-2016-fv6.0.nc';
ncfile_2015 = 'ESACCI-OC-L3S-CHLOR_A-MERGED-1Y_YEARLY_4km_GEO_PML_OCx-2015-fv6.0.nc';    
ncfile_2014 = 'ESACCI-OC-L3S-CHLOR_A-MERGED-1Y_YEARLY_4km_GEO_PML_OCx-2014-fv6.0.nc';
ncfile_2013 = 'ESACCI-OC-L3S-CHLOR_A-MERGED-1Y_YEARLY_4km_GEO_PML_OCx-2013-fv6.0.nc';
ncfile_2012 = 'ESACCI-OC-L3S-CHLOR_A-MERGED-1Y_YEARLY_4km_GEO_PML_OCx-2012-fv6.0.nc';
ncfile_2011 = 'ESACCI-OC-L3S-CHLOR_A-MERGED-1Y_YEARLY_4km_GEO_PML_OCx-2011-fv6.0.nc';
ncfile_2010 = 'ESACCI-OC-L3S-CHLOR_A-MERGED-1Y_YEARLY_4km_GEO_PML_OCx-2010-fv6.0.nc';
ncfile_2009 = 'ESACCI-OC-L3S-CHLOR_A-MERGED-1Y_YEARLY_4km_GEO_PML_OCx-2009-fv6.0.nc';
ncfile_2008 = 'ESACCI-OC-L3S-CHLOR_A-MERGED-1Y_YEARLY_4km_GEO_PML_OCx-2008-fv6.0.nc';
ncfile_2007 = 'ESACCI-OC-L3S-CHLOR_A-MERGED-1Y_YEARLY_4km_GEO_PML_OCx-2007-fv6.0.nc';
ncfile_2006 = 'ESACCI-OC-L3S-CHLOR_A-MERGED-1Y_YEARLY_4km_GEO_PML_OCx-2006-fv6.0.nc';
ncfile_2005 = 'ESACCI-OC-L3S-CHLOR_A-MERGED-1Y_YEARLY_4km_GEO_PML_OCx-2005-fv6.0.nc';
ncfile_2004 = 'ESACCI-OC-L3S-CHLOR_A-MERGED-1Y_YEARLY_4km_GEO_PML_OCx-2004-fv6.0.nc';
ncfile_2003 = 'ESACCI-OC-L3S-CHLOR_A-MERGED-1Y_YEARLY_4km_GEO_PML_OCx-2003-fv6.0.nc';
ncfile_2002 = 'ESACCI-OC-L3S-CHLOR_A-MERGED-1Y_YEARLY_4km_GEO_PML_OCx-2002-fv6.0.nc';

% Load lat, lon, and chlor_a data
lat = ncread(ncfile_2022, 'lat');
lon = ncread(ncfile_2022, 'lon');
chl_2022 = ncread(ncfile_2022, 'chlor_a');
chl_2021 = ncread(ncfile_2021, 'chlor_a');
chl_2020 = ncread(ncfile_2020, 'chlor_a');
chl_2019 = ncread(ncfile_2019, 'chlor_a');
chl_2018 = ncread(ncfile_2018, 'chlor_a');
chl_2017 = ncread(ncfile_2017, 'chlor_a');
chl_2016 = ncread(ncfile_2016, 'chlor_a');
chl_2015 = ncread(ncfile_2015, 'chlor_a');
chl_2014 = ncread(ncfile_2014, 'chlor_a');
chl_2013 = ncread(ncfile_2013, 'chlor_a');
chl_2012 = ncread(ncfile_2012, 'chlor_a');
chl_2011 = ncread(ncfile_2011, 'chlor_a');
chl_2010 = ncread(ncfile_2010, 'chlor_a');
chl_2009 = ncread(ncfile_2009, 'chlor_a');
chl_2008 = ncread(ncfile_2008, 'chlor_a');
chl_2007 = ncread(ncfile_2007, 'chlor_a');
chl_2006 = ncread(ncfile_2006, 'chlor_a');
chl_2005 = ncread(ncfile_2005, 'chlor_a');
chl_2004 = ncread(ncfile_2004, 'chlor_a');
chl_2003 = ncread(ncfile_2003, 'chlor_a');
chl_2002 = ncread(ncfile_2002, 'chlor_a');

%switch lonxlat to latxlon
chl_2022 = chl_2022';
chl_2021 = chl_2021';
chl_2020 = chl_2020';
chl_2019 = chl_2019';
chl_2018 = chl_2018';
chl_2017 = chl_2017';
chl_2016 = chl_2016';
chl_2015 = chl_2015';
chl_2014 = chl_2014';
chl_2013 = chl_2013';
chl_2012 = chl_2012';
chl_2011 = chl_2011';
chl_2010 = chl_2010';
chl_2009 = chl_2009';
chl_2008 = chl_2008';
chl_2007 = chl_2007';
chl_2006 = chl_2006';
chl_2005 = chl_2005';
chl_2004 = chl_2004';
chl_2003 = chl_2003';
chl_2002 = chl_2002';

% Set up map projection
figure;
m_proj('miller','lon',[-180 180],'lat',[-60 60]);
% Add 2021 annual mean chlorophyll
m_pcolor(glob_lon, glob_lat, chl)
% Add chl=0.1 contours
hold on
h1 = m_contour(lon, lat, chl_2022,[0.1 0.1],'color',[0.5 0.5 0.5],'LineWidth',0.5);
m_contour(lon, lat, chl_2021,[0.1 0.1],'color',[0.5 0.5 0.5],'LineWidth',0.5);
m_contour(lon, lat, chl_2020,[0.1 0.1],'color',[0.5 0.5 0.5],'LineWidth',0.5);
m_contour(lon, lat, chl_2019,[0.1 0.1],'color',[0.5 0.5 0.5],'LineWidth',0.5);
m_contour(lon, lat, chl_2018,[0.1 0.1],'color',[0.5 0.5 0.5],'LineWidth',0.5);
m_contour(lon, lat, chl_2017,[0.1 0.1],'color',[0.5 0.5 0.5],'LineWidth',0.5);
m_contour(lon, lat, chl_2016,[0.1 0.1],'color',[0.5 0.5 0.5],'LineWidth',0.5);
m_contour(lon, lat, chl_2015,[0.1 0.1],'color',[0.5 0.5 0.5],'LineWidth',0.5);
m_contour(lon, lat, chl_2014,[0.1 0.1],'color',[0.5 0.5 0.5],'LineWidth',0.5);
m_contour(lon, lat, chl_2013,[0.1 0.1],'color',[0.5 0.5 0.5],'LineWidth',0.5);
m_contour(lon, lat, chl_2012,[0.1 0.1],'color',[0.5 0.5 0.5],'LineWidth',0.5);
m_contour(lon, lat, chl_2011,[0.1 0.1],'color',[0.5 0.5 0.5],'LineWidth',0.5);
m_contour(lon, lat, chl_2010,[0.1 0.1],'color',[0.5 0.5 0.5],'LineWidth',0.5);
m_contour(lon, lat, chl_2009,[0.1 0.1],'color',[0.5 0.5 0.5],'LineWidth',0.5);
m_contour(lon, lat, chl_2008,[0.1 0.1],'color',[0.5 0.5 0.5],'LineWidth',0.5);
m_contour(lon, lat, chl_2007,[0.1 0.1],'color',[0.5 0.5 0.5],'LineWidth',0.5);
m_contour(lon, lat, chl_2006,[0.1 0.1],'color',[0.5 0.5 0.5],'LineWidth',0.5);
m_contour(lon, lat, chl_2005,[0.1 0.1],'color',[0.5 0.5 0.5],'LineWidth',0.5);
m_contour(lon, lat, chl_2004,[0.1 0.1],'color',[0.5 0.5 0.5],'LineWidth',0.5);
m_contour(lon, lat, chl_2003,[0.1 0.1],'color',[0.5 0.5 0.5],'LineWidth',0.5);
m_contour(lon, lat, chl_2003,[0.1 0.1],'color',[0.5 0.5 0.5],'LineWidth',0.5);
m_contour(glob_lon, glob_lat, annual_mean(:,:,years==1997),[0.1 0.1],'color',[0.5 0.5 0.5],'LineWidth',0.5)
m_contour(glob_lon, glob_lat, annual_mean(:,:,years==1998),[0.1 0.1],'color',[0.5 0.5 0.5],'LineWidth',0.5)
m_contour(glob_lon, glob_lat, annual_mean(:,:,years==1999),[0.1 0.1],'color',[0.5 0.5 0.5],'LineWidth',0.5)
m_contour(glob_lon, glob_lat, annual_mean(:,:,years==2000),[0.1 0.1],'color',[0.5 0.5 0.5],'LineWidth',0.5)
m_contour(glob_lon, glob_lat, annual_mean(:,:,years==2001),[0.1 0.1],'color',[0.5 0.5 0.5],'LineWidth',0.5)
% Add hand-drawn masks I use...
h2 = m_line(NA_lon_chl, NA_lat_chl, 'color','k','linewi',2)
m_line(SA_lon_chl, SA_lat_chl, 'color','k','linewi',2)
m_line(I_lon_chl, I_lat_chl, 'color','k','linewi',2)
m_line(NP1_lon_chl, NP1_lat_chl, 'color','k','linewi',2)
m_line(SP1_lon_chl, SP1_lat_chl, 'color','k','linewi',2)
m_line(NP2_lon_chl, NP2_lat_chl, 'color','k','linewi',2)
m_line(SP2_lon_chl, SP2_lat_chl, 'color','k','linewi',2)
c = colorbar
cmocean('algae')
caxis([0 0.3])
c.Label.String = 'Chl-a concentration (mg m^-^3)'G
m_grid('fontsize',9);
m_coast('color','k')
set(gcf,'color','w')
set(gca,'Fontsize',12)
hold off;

legend([h1, h2], {'1997-2022', '2021'}, 'Location', 'south')