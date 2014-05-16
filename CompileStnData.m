
clear;
clc;

if strcmpi(getenv('OS'),'Windows_NT')
    addpath('C:\Users\gianotti\Documents\IntensityLib\')
else % linux
    addpath('IntensityLib');
end

today_doy = today - datenum(year(today),1,1) + 1;
remove_seasonal_cycle = false;

%% Download and calculate all of the LL data:
calculate_daily_LL_data;

%% Create different accumulated LL time series (each gets saved with the 
% daily LL data for that stn).

% Delete old accumulated data files:
delete('*accum_DOY*.mat');
accumulate_LL({'oneyear','twoyear','threeyear'},today_doy,remove_seasonal_cycle);
accumulate_LL({'oneyear','twoyear','threeyear'},1,remove_seasonal_cycle);

% Transform the annualized (or bi-annualized, etc.) data to be more normal
% using the sim distribution
cdf_transform_accum_LL_data();

% Cluster the obs data and make some plots:
make_accumlated_cluster_plots_LL(1);
make_accumlated_cluster_plots_LL(today_doy);

%% Now do the same thing for total precipitation and occurrence instead of LL:

accumulate_precip({'oneyear','twoyear','threeyear'},today_doy,remove_seasonal_cycle);
accumulate_precip({'oneyear','twoyear','threeyear'},1,remove_seasonal_cycle);

% Transform the annualized (or bi-annualized, etc.) data to be more normal 
% using the sim distribution
cdf_transform_accum_precip_data();

% Cluster the obs data and make some plots:
make_accumlated_cluster_plots(1);
make_accumlated_cluster_plots(today_doy);

% Make annual LL/precip scatter plots and seasonal time-series:
make_annual_LL_precip_plots;

%% Try a comparison between May 1 start date with and without May-Oct precip:
may1 = 121;
oct31 = 304;
accumulate_LL({'oneyear','twoyear','threeyear'},may1,0);
accumulate_precip({'oneyear','twoyear','threeyear'},may1,0);

accumulate_precip_wet_season_only({'oneyear','twoyear','threeyear'},may1,0);
accumulate_LL_wet_season_only({'oneyear','twoyear','threeyear'},may1,0);

cdf_transform_accum_LL_data;
cdf_transform_accum_precip_data;

make_accumulated_cluster_plots_wet_season(may1);


