function [LL, occ_frac] = get_log_likelihood89(stn_id, simobs, season)
% The function get_log_likelihood89 returns the log likelihood of the 
% 89-day precipitation record for a given station, in a particular 89-day 
% season, using either the observed or simulated record. 
% 
% The input stn_id is a string of the 6-digit USHCN id for the station.
%
% The input simobs is a string, either 'sim' or 'obs', determining which
% data to use in calculating the Log Likelihood. If simobs is 'obs', then
% the output of this function will be a n_year x 1 vector of log 
% likelihoods. If simobs is 'sim', then the output is an n_year x 1000 
% matrix of log likelihoods.
%
% The input season is one of 'DJF', 'MAM', 'JJA', or 'SON'.
% 
% The second output, occ_frac, is a number between 0 and 1 that says what
% fraction of the LL is due to occurrence versus intensity.


switch season
    case {'DJF','djf'}
        center_date = 15;
    case {'MAM','mam'}
        center_date = 105;
    case {'JJA','JJA'}
        center_date = 196;
    case {'SON','SON'}
        center_date = 288;    
end

if strcmp(simobs,'obs')
    % obs
    stn = load_stn_data(stn_id,'ImpStn');
    % Shift the data so that we only have the season in question and the
    % previous 5 days (for occ state):
    data = ShiftXdays(stn.intensity_data, 50-center_date);
    data = data(2:(end-1), 1:94); % 89 days, plus 5 in the begining to determine occ state
    
    n_years = size(data,1); % going to remove first and last year
else
    % sim
    stn = load_stn_data(stn_id,'SimStn');
    % Shift the data so that we only have the season in question and the
    % previous 5 days (for occ state):
    data = ShiftXdays(stn.intensity_data, 50-center_date);
    data = data(:, 1:94); % 89 days, plus 5 in the begining to determine occ state

    n_years = size(data,1);
end

LL_occ = nan(n_years,89);
LL_int = nan(n_years,89);

% Load the int model data:
int_model_path = ['C:\Users\gianotti\Documents\IntensityData\IntensityModelData\IntensityModelData',stn_id,'.mat'];
load(int_model_path);

% load the occ model data:
occ_mod = load_stn_data(stn_id,'SelectedAICMod');

for day = 1:89
    day_of_year = mod(center_date-46+day,365)+1;
    occ_history = (data(:, day:(day+4)) > 0); % an n_years x 5 mtx of 0s and 1s for occurrence history
    % Going to do this diredctly from history. Don't need to translate
    % twice...
    %occ_states = get_occ_states(occ_history); % an n_years x 1 vector of occ states, between 1 and 32
    rain_today = (data(:,day+5) > 0); % n_years x 1 vector of occurrence on THIS day
    % Determine probability of occurrence (given history)
    p_occ_today = get_prob_occurrence(occ_history,occ_mod.trans_probs{day_of_year});
    % Assign LL for occurrence
    LL_occ(:,day) = log(( (rain_today == 1) .* p_occ_today) + ((rain_today == 0) .* (1-p_occ_today)));
    % Assign LL for intensity
    LL_int(:,day) = assign_LL_int(data(:,day+5),best_params{day});
    
end

% You can add up the days first, or add the int and occ first, it doesn't
% matter. We're going to add up the days first so that we can attribute the
% LL to occ and int appropriately.

% Now add up LL for all days:
LL_occ = sum(LL_occ,2);
LL_int = sum(LL_int,2);

LL = LL_occ + LL_int;

occ_frac = LL_occ ./ LL;

if strcmp(simobs,'sim')
    LL = reshape(LL,[numel(LL)/1000,1000]);
    occ_frac = reshape(occ_frac, [numel(LL)/1000,1000]);
end

end % function