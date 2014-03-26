function [LL, years, occ_frac] = get_daily_log_likelihood(stn_id, data, start_year)
% The function get_daily_log_likelihood determines the log likelihood of
% the data (input 2) given the Climate-Stationary Weather Model (CSWM) for
% station id. Inputs are:
%
% stn_id: the 6-character GHCN id string,
%
% data: an N x 365 matrix of precipitation data INCLUDING missing values as
% NaN,
%
% start_year: the (integer) year corresponding to the first year of
% observed data,
%
% normalize_LL: a binary (true/false) flag for whether to divide the LL
% data by the number of observed data points for that year.
%
% Outputs are:
%
% LL: an N x 365 vector of log-likelihoods, and
%
% years: an N x 1 vector of the associated years
%
% occ_frac: the fraction of the log-likelihood that is due to intensity

% Check to make sure data is N x 365:
if (size(data,2) ~= 365)
    error('The input data to the function get_annual_log_likelihood should be an N x 365 matrix. Aborting!');
end

n_years = size(data,1);

LL_occ = nan(n_years,365);
LL_int = nan(n_years,365);

% Load the int model data:
int_model_path = ['C:\Users\gianotti\Documents\IntensityData\IntensityModelData\IntensityModelData',stn_id,'.mat'];
load(int_model_path);

% load the occ model data:
occ_mod = load_stn_data(stn_id,'SelectedAICMod');


for day = 1:365
    day
    % Shift the data so that 5-day history is columns 1-5, today is column
    % 6 (we'll need to throw out year 1):
    shifted_data = ShiftXdays(data, 6 - day);
    history = shifted_data(:,1:5);
    todays_rain = shifted_data(:,6);
  
    % Probability of occurrence today
    p_occ_today = get_prob_occurrence( (history>0), occ_mod.trans_probs{day});
  
    % Assign LL for occurrence
    LL_occ(:,day) = log(( (todays_rain > 0) .* p_occ_today) + ((todays_rain == 0) .* (1-p_occ_today)));
    
    % Assign LL for intensity     
    int_chain_order = log2(size(best_params{day},2));
    switch int_chain_order
        case 0
            LL_int(:,day) = assign_LL_int(todays_rain, best_params{day});
        case 1 
            params = best_params{day};
            rained_yesterday = (history(:,5) > 1);
            no_rain_yesterday = (history(:,5) == 0); % Need seperate cases due to nans
            nan_yesterday = isnan(history(:,5));
            LL_int(rained_yesterday,day) = assign_LL_int(todays_rain(rained_yesterday), params(:,2));
            LL_int(no_rain_yesterday,day) = assign_LL_int(todays_rain(no_rain_yesterday), params(:,1));
            LL_int(nan_yesterday,day) = nan;

        otherwise
            error('We need to add cases for higher chain orders still. Aborting!');
    end   
end

% You can add up the days first, or add the int and occ first, it doesn't
% matter. We're going to add up the days first so that we can attribute the
% LL to occ and int appropriately.

% Now add up LL for all days (but exclude the first year):
LL_occ(1,:) = []; % Throw out first year...
LL_int(1,:) = [];

% throw out bad data:
bad_data = (isnan(LL_occ) | isnan(LL_int));
LL_occ(bad_data) = nan;
LL_int(bad_data) = nan;

LL = LL_occ + LL_int;

occ_frac = LL_occ ./ LL;


years = ((start_year+1):(start_year+n_years-1))';
end % function
