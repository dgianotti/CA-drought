% Let's get the (current) data for the CA stations!

clear;
clc;

addpath('C:\Users\gianotti\Documents\IntensityLib');

load('CA_ids.mat');

fraction_missing_vec = [ 0.2302, 0.3140, 1, 0.2133, 0.4727, ...
    0.0111, 0.0364, 0.0078, 0.2055, 0.0559, ...
    0.6050, 0.0013, 0.0124, 0.3296, 0.2581, ...
    0.1248, 0.0377, 0.0013, 0.2302, 0.0663, ...
    0.0026, 0.0020, 0.0104, 0.0228, 0.0449, ...
    0.1990, 0.5379, 0.1586, 0.1944, 0.5879, ...
    1, 0.0845, 0.1118, 0.1190, 0.0059 ];

good_CA_IDs = CA_IDs(fraction_missing_vec < 0.05);
    

%% Loop over stations:
for i = 1:length(good_CA_IDs)

    id = good_CA_IDs{i}
    % Download the latest USHCN data:
    url = ['http://www1.ncdc.noaa.gov/pub/data/ghcn/daily/hcn/USC00',id,'.dly'];
    filename = ['GHCN-Daily/',id,'.dly'];
    %urlwrite(url,filename);
    
    
    % Read in fixed-width data:
    [datenums,precip] = format_GHCN_precip_data(filename,'PadFirstLastYears','ExcludeLeapDays');

    %     % Set last day to to:
    %     feb28_2014 = datenum('2014-02-28','yyyy-mm-dd');
    %     datenums = [datenums(:); ((datenums(end)+1):feb28_2014)'];
    %     precip = [precip(:); nan([length(datenums)-length(precip),1])];
        
    new_precip = precip(datenums>= datenum('2010-01-01') & datenums<= today);
    new_datenums = datenums(datenums>= datenum('2010-01-01') & datenums<= today);

    % Determine fraction of data missing since 2010:
    fraction_missing = sum(isnan(new_precip))/numel(new_precip)

    
    %if (mod(length(new_precip),365) ~= 0) %uh oh... it should...
    %error('The new data should be a multiple of 365 days... did you use ExcludeLeapDays in format_GHCN_precip_data\n');
    %end


    % Fill missing data values somehow!
    % Neighboring GHCN sites?
    missing_data_datenums = new_datenums(isnan(new_precip)) ;
    
    % Look up lat/lon in ushcn-stations.txt:
    [lat, lon] = get_ushcn_lat_lon(id);
    
    % Maybe use GSOD:
    % ftp://ftp.ncdc.noaa.gov/pub/data/gsod
    
    % filled_precip = get_precip_GSOD(missing_data_datenums,lat,lon,20);

    % Or NOAA's Quality controlled local climatological data:
    % http://cdo.ncdc.noaa.gov/qclcd/QCLCD
    % filled_precip = get_precip_QCLCD(missing_data_datenums,lat,lon);

    % Load old precip data:
    ImpStn = load_stn_data(id,'ImpStn');

    % Determine likelihood of data given model:
    
    nan_padded_new_data = nan(365*5,1);
    nan_padded_new_data(1:length(new_precip)) = new_precip;
    
    data = [ImpStn.intensity_data; reshape(nan_padded_new_data,[365,5])'];
        
    normalize_LL = true;
    start_year = 2010 - ImpStn.num_years;
    [LL_obs, years] = get_annual_log_likelihood(id, data, start_year, normalize_LL);
    
    % Now sim data!
    SimStn = load_stn_data(id,'SimStn');        

    [LL_sim, ~] = get_annual_log_likelihood(id, SimStn.intensity_data, start_year, normalize_LL);

    n_years = length(years);
    num_useable_sims = floor(length(LL_sim)/n_years);
    LL_sim = LL_sim(1:(n_years*num_useable_sims));
    LL_sim = reshape(LL_sim, [n_years,num_useable_sims]);
    
    plot(years,LL_obs,'-k');
    hold on;
    plot(years,quantile(LL_sim,.05,2),'-b');
    plot(years,quantile(LL_sim,.95,2),'-b');
    plot(years,quantile(LL_sim,.5,2),'-r');
    ylabel('Log Likelihood');
    title(sprintf('Normalized annual log likelihood precipitation for station %s',id)); 
    
    print(gcf,'-dpng',sprintf('LL_%s.png',id));
    save(sprintf('LL_%s.mat',id),'LL_obs','LL_sim','years');
    
end 
    
    
    