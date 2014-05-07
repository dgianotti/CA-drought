function accumulate_LL(periods, start_DOY, remove_seasonal_cycle)
% This will save the accumulated LL time series in the already created
% LL_[id].mat files.
%
% periods should be a cell array of strings (allowable strings so far are:
%
% 'oneyear','twoyear','threeyear'
% start_DOY should be an integer between 1 and 365.

load('CA_ids.mat');

% Loop over each station
for i = 1:length(good_CA_IDs)
    fprintf('Accumulating likelihoods for station %i of %i...\n',i,length(good_CA_IDs));
    id = good_CA_IDs{i};

    % Load the daily LL data:
    filename = ['LL_',id,'.mat'];
    load(filename);
    
    % Remove the seasonal cycle if desired:
    if remove_seasonal_cycle
        LL_mean = nanmean(LL_sim,1);
        LL_std = nanstd(LL_sim,0,1);
        LL_obs = (LL_obs - repmat(LL_mean, [size(LL_obs,1),1]))...
            ./ repmat(LL_std, [size(LL_obs,1),1]);
        LL_sim = (LL_sim - repmat(LL_mean, [size(LL_sim,1),1]))...
            ./ repmat(LL_std, [size(LL_sim,1),1]);
    end
    
    % Shift the data so that it starts on the given start_DOY:
    LL_obs = ShiftXdays(LL_obs,1-start_DOY);
    LL_obs(end,:) = []; % This one has both the incomplete first year and the NaNs from the end of this year.
    LL_sim = ShiftXdays(LL_sim,1-start_DOY);
    LL_sim(end,:) = [];
    years(1) = [];

    % % % % % % % % % % % % % % %
    % 1yr
    if any(strcmpi('oneyear',periods)) || any(strcmpi('1yr',periods))

        LL_obs_tmp = LL_obs;
        LL_sim_tmp = LL_sim;
        
        % Now we need to make the sim data an integer multiple of the obs data
        % and insert NaNs at the correct places:
        n_years = size(LL_obs_tmp,1);
        n_sims = floor(size(LL_sim_tmp,1)/n_years);
        LL_sim_tmp( (n_years*n_sims+1):end, : ) = [];
        obs_nans = isnan(LL_obs_tmp);
        LL_sim_tmp( repmat(obs_nans,[n_sims,1]) ) = nan;
        
        LL_obs_1yr = nanmean(LL_obs_tmp,2); % n_years x 1
        LL_sim_1yr = reshape(nanmean(LL_sim_tmp,2), [n_years,n_sims]); % As an n_years x n_sims matrix        
        
    end % if
    
    % % % % % % % % % % % % % % %
    % 2yr
    if any(strcmpi('twoyear',periods)) || any(strcmpi('2yr',periods))
        LL_obs_tmp = LL_obs;
        LL_sim_tmp = LL_sim;
        
        % Make the LL matrices have an even number of years:
        years_to_remove = mod(size(LL_obs,1),2);
        LL_obs_tmp(1:years_to_remove,:) = []; % Even!
        
        % Now we need to make the sim data an integer multiple of the obs data
        % and insert NaNs at the correct places:
        n_years = size(LL_obs_tmp,1);
        n_sims = floor(size(LL_sim_tmp,1)/n_years);
        LL_sim_tmp( (n_years*n_sims+1):end, : ) = [];
        obs_nans = isnan(LL_obs_tmp);
        LL_sim_tmp( repmat(obs_nans,[n_sims,1]) ) = nan;

        % Now make each row represent 2 years:
        LL_obs_tmp = reshape(LL_obs_tmp',[2*365, numel(LL_obs_tmp)/(2*365)])';
        LL_sim_tmp = reshape(LL_sim_tmp',[2*365, numel(LL_sim_tmp)/(2*365)])';
        
        LL_obs_2yr = nanmean(LL_obs_tmp,2); % n_years x 1
        LL_sim_2yr = reshape(nanmean(LL_sim_tmp,2),... 
            [length(LL_obs_2yr),size(LL_sim_tmp,1)/length(LL_obs_2yr)]); % As an n_years x n_sims matrix        
    end % if for 2yr
    
    
    % % % % % % % % % % % % % % %
    % 3yr
    if any(strcmpi('threeyear',periods)) || any(strcmpi('3yr',periods))
        LL_obs_tmp = LL_obs;
        LL_sim_tmp = LL_sim;
        
        % Make the LL matrices have a number of years that si a multiple of 3:
        years_to_remove = mod(size(LL_obs,1),3);
        LL_obs_tmp(1:years_to_remove,:) = []; % Multiple of 3!
        
        % Now we need to make the sim data an integer multiple of the obs data
        % and insert NaNs at the correct places:
        n_years = size(LL_obs_tmp,1);
        n_sims = floor(size(LL_sim_tmp,1)/n_years);
        LL_sim_tmp( (n_years*n_sims+1):end, : ) = [];
        obs_nans = isnan(LL_obs_tmp);
        LL_sim_tmp( repmat(obs_nans,[n_sims,1]) ) = nan;

        % Now make each row represent 3 years:
        LL_obs_tmp = reshape(LL_obs_tmp',[3*365, numel(LL_obs_tmp)/(3*365)])';
        LL_sim_tmp = reshape(LL_sim_tmp',[3*365, numel(LL_sim_tmp)/(3*365)])';
        
        LL_obs_3yr = nanmean(LL_obs_tmp,2); % n_years x 1
        LL_sim_3yr = reshape(nanmean(LL_sim_tmp,2),... 
            [length(LL_obs_3yr),size(LL_sim_tmp,1)/length(LL_obs_3yr)]); % As an n_years x n_sims matrix        
    end % if for 3yr

    % Now save the new LL time series:
    % For now, let's just assume that they've all been created:
    filename = sprintf('LL_%s_accum_DOY%d.mat',id,start_DOY);
    save(filename, 'LL_obs_1yr','LL_sim_1yr', 'LL_obs_2yr','LL_sim_2yr',... 
        'LL_obs_3yr','LL_sim_3yr','id','start_DOY');
    
end




end