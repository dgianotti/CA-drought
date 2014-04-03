function calculate_daily_LL_data()
% This function will download/update all of the USHCN data for our selected
% stations, and recalculate the daily log-likelihoods. 
%
% Each file, saved as LL_id.mat, will contain an N x 1 vector, 'years', an
% N x 365 matrix 'LL_obs', another the same size called 'LL_obs_std_norm'
% (which has undergone a preliminary attempt at normalization, but still
% has a seasonal cycle), a roughly-but-not-exactly (996*N) x 365 matrix
% 'LL_sim', and another, 'LL_sim_std_norm'.
%
% All data will still have a seasonal likelihood cycle (due to the seasonal
% cycle of occurrence and its effects on likelihood).


%addpath('C:\Users\gianotti\Documents\IntensityLib');

load('CA_ids.mat');

fraction_missing_vec = [ 0.2302, 0.3140, 1, 0.2133, 0.4727, ...
    0.0111, 0.0364, 0.0078, 0.2055, 0.0559, ...
    0.6050, 0.0013, 0.0124, 0.3296, 0.2581, ...
    0.1248, 0.0377, 0.0013, 0.2302, 0.0663, ...
    0.0026, 0.0020, 0.0104, 0.0228, 0.0449, ...
    0.1990, 0.5379, 0.1586, 0.1944, 0.5879, ...
    1, 0.0845, 0.1118, 0.1190, 0.0059 ];

good_CA_IDs = CA_IDs(fraction_missing_vec < 0.05);

% save('CA_ids.mat','CA_IDs','good_CA_IDs');

% Loop over stations:
for i = 1:length(good_CA_IDs)

    id = good_CA_IDs{i};
    fprintf('Downloading data for station %s...\n',id);
    % Download the latest USHCN data:
    url = ['http://www1.ncdc.noaa.gov/pub/data/ghcn/daily/hcn/USC00',id,'.dly'];
    filename = ['GHCN-Daily/',id,'.dly'];
    urlwrite(url,filename);
    
    
    % Read in fixed-width data:
    fprintf('Re-formatting data for station %s...\n',id);
    [datenums,precip] = format_GHCN_precip_data(filename,'PadFirstLastYears','ExcludeLeapDays');

    %     % Set last day to to:
    %     feb28_2014 = datenum('2014-02-28','yyyy-mm-dd');
    %     datenums = [datenums(:); ((datenums(end)+1):feb28_2014)'];
    %     precip = [precip(:); nan([length(datenums)-length(precip),1])];
        
    new_precip = precip(datenums>= datenum('2010-01-01') & datenums<= today);
    new_datenums = datenums(datenums>= datenum('2010-01-01') & datenums<= today);

    % Determine fraction of data missing since 2010:
    %fraction_missing = sum(isnan(new_precip))/numel(new_precip);
    
    %missing_data_datenums = new_datenums(isnan(new_precip));
    
    % Look up lat/lon in ushcn-stations.txt:
    % [lat, lon] = get_ushcn_lat_lon(id);
    
    % Load old precip data:
    ImpStn = load_stn_data(id,'ImpStn');

    % Determine likelihood of data given model:
    
    nan_padded_new_data = nan(365*5,1);
    nan_padded_new_data(1:length(new_precip)) = new_precip;
    
    data = [ImpStn.intensity_data; reshape(nan_padded_new_data,[365,5])'];
        
    start_year = 2010 - ImpStn.num_years;
    fprintf('Calculating daily observed LL for station %s...\n',id);
    [LL_obs_std_norm,~] = get_daily_log_likelihood_std_normal(id, data, start_year);
    [LL_obs, years] = get_daily_log_likelihood(id, data, start_year);
    
    % Now sim data!
    SimStn = load_stn_data(id,'SimStn');        

    fprintf('Calculating daily simulated LL for station %s...\n',id);
    [LL_sim_std_norm, ~] = get_daily_log_likelihood_std_normal(id, SimStn.intensity_data, start_year);
    [LL_sim, ~] = get_daily_log_likelihood(id, SimStn.intensity_data, start_year);
    save(sprintf('LL_%s.mat',id),'LL_obs','LL_sim','LL_obs_std_norm','LL_sim_std_norm','years');
    
end 