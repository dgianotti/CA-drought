% Let's get the (current) data for the CA stations!

clear;
clc;

%addpath('C:\Users\gianotti\Documents\IntensityLib');

load('CA_ids.mat');

%% Download data:
% for i = 1:length(CA_IDs)
id = CA_IDs{1}; % Berkeley
% Download the latest USHCN data:
url = ['http://www1.ncdc.noaa.gov/pub/data/ghcn/daily/hcn/USC00',id,'.dly'];
filename = [id,'.dly'];
urlwrite(url,filename);

%% Read in fixed-width data:

[datenums,precip] = format_GHCN_precip_data(filename,'PadFirstLastYears','ExcludeLeapDays','MFLAG_P_AsNaN');

new_precip = precip(datenums>= datenum('2010-01-01'));
new_datenums = datenums(datenums>= datenum('2010-01-01'));

if mod(length(new_precip),365) ~= 0 %uh oh... it should...
    error('The new data should be a multiple of 365 days... did you use ExcludeLeapDays in format_GHCN_precip_data\n');
end


%% Fill missing data values somehow!
% Neighboring GHCN sites for Berkeley:

% Maybe use GSOD:
% ftp://ftp.ncdc.noaa.gov/pub/data/gsod
missing_data_datenums = new_datenums(isnan(new_precip)) ;

lat = 37 + (52/60) + (18/3600);
lon = -(122 + (16/60) + (22/3600));

filled_precip = get_precip_GSOD(missing_data_datenums,lat,lon,7.5);

% Or NOAA's Quality controlled local climatological data:
% http://cdo.ncdc.noaa.gov/qclcd/QCLCD

%% Load old precip data:
ImpStn = load_stn_data(id,'ImpStn');

intensity_data = ImpStn.intensity_data;


new_intensity_data = [intensity_data; reshape(new_precip,[length(new_precip)/365,365])];
sum(isnan(new_precip))

%% Get model data


%% Determine likelihood of data given model:

%% Determine