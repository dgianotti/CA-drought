function [LL_std_norm, years] = get_daily_log_likelihood_std_normal(stn_id, data, start_year)
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

% Load the int model data:
int_model_path = ['C:\Users\gianotti\Documents\IntensityData\IntensityModelData\IntensityModelData',stn_id,'.mat'];
load(int_model_path);

% load the occ model data:
occ_mod = load_stn_data(stn_id,'SelectedAICMod');

addpath('C:\Users\gianotti\Documents\IntensityLib');

LL_std_norm = zeros(size(data,1),365);

for day = 1:365
    fprintf('%i\n',day);

    occ_chain_order = occ_mod.selected_orders(day);
    int_chain_order = log2(size(best_params{day},2));
    
    % Shift the data so that 5-day history is columns 1-5, today is column
    % 6 (we'll need to throw out year 1):
    shifted_data = ShiftXdays(data, 6 - day);

    missing_data = isnan(shifted_data);
    % Set missing data = 0, then we need to fix it at the end
    shifted_data(missing_data) = 0;

    history = shifted_data(:,1:5);
    todays_rain = shifted_data(:,6);
    
    % Probability of occurrence today
    p_occ_today = get_prob_occurrence( (history>0), occ_mod.trans_probs{day});
  
    cdfs = zeros(size(todays_rain));
       
    params = best_params{day};
    
    % Assign LL for intensity     
    switch int_chain_order
        case 0
            cdfs(todays_rain == 0) = 1 - p_occ_today(todays_rain == 0);
            cdfs(todays_rain > 0) = 1 - p_occ_today(todays_rain > 0)...
                + p_occ_today(todays_rain > 0).* mixgamcdf(todays_rain(todays_rain > 0),params);                        
        case 1 
            rained_yesterday = (history(:,5) > 1);
            
            % Four cases: rained yesterday and today (11), yesterday but
            % not today (10), not yesterday but today (01), neither (00):
            case_11 = rained_yesterday & (todays_rain > 0); % use params(:,2)
            case_10 = rained_yesterday & (todays_rain == 0); 
            case_01 = ~rained_yesterday & (todays_rain > 0); % use params(:,1)
            case_00 = ~rained_yesterday & (todays_rain == 0); 
            
            cdfs(case_11) = 1 - p_occ_today(case_11)...
                + p_occ_today(case_11).* mixgamcdf(todays_rain(case_11),params(:,2));
            cdfs(case_10) = 1 - p_occ_today(case_10);
            cdfs(case_01) = 1 - p_occ_today(case_01)...
                + p_occ_today(case_01).* mixgamcdf(todays_rain(case_01),params(:,1));
            cdfs(case_00) = 1 - p_occ_today(case_00); 
            
        otherwise
            error('We need to add cases for higher chain orders still. Aborting!');
    end   
                   
    % Transform to normal using inverse normal CDF:
    LL_std_norm(:,day) = log(normpdf(norminv(cdfs,0,1),0,1)); % norminv(cdfs) is a precip equivalent, normlike(norminv(cdfs)) is -LL
    
    
    % Now need to fix missing data:
    max_chain_order = max(occ_chain_order, int_chain_order);
    for j = 0:max_chain_order
        LL_std_norm( missing_data(:,6-j) ) = nan;
    end
    
end

% remove the first year:
LL_std_norm(1,:) = [];


years = ((start_year+1):(start_year+n_years-1))';
end % function
