function accumulate_precip(periods, start_DOY, remove_seasonal_cycle)
% This will save the accumulated precip time series in precip_[id].mat 
% files.
%
% periods should be a cell array of strings (allowable strings so far are:
%
% 'oneyear','twoyear','threeyear'
% start_DOY should be an integer between 1 and 365.

load('CA_ids.mat');
load('nino34.mat');
nino34_monthly(nino34_monthly == -99.99) = nan;
nino_djfma = ShiftXdays(nino34_monthly,1);
nino_djfma = nanmean(nino_djfma(:,1:5),2);
nino_djfma(end) = []; % Throw out 2014 (IS THIS CORRECT TO DO?!?!?!?!)		

% Loop over each station
for i = 1:length(good_CA_IDs)
    fprintf('Accumulating precipitation for station %i of %i...\n',i,length(good_CA_IDs));
    id = good_CA_IDs{i};

    % Load the daily LL data:
    filename = ['LL_',id,'.mat'];
    load(filename);
    
    precip_obs(1,:) = [];
    
    % Load the sim data too:
    SimStn = load_stn_data(id,'SimStn');
    precip_sim = SimStn.intensity_data;
    
    
    % Remove the seasonal cycle if desired:
    if remove_seasonal_cycle
        precip_mean = nanmean(precip_sim,1);
        precip_std = nanstd(precip_sim,0,1);
        precip_obs = (precip_obs - repmat(precip_mean, [size(precip_obs,1),1]))...
            ./ repmat(precip_std, [size(precip_obs,1),1]);
        precip_sim = (precip_sim - repmat(precip_mean, [size(precip_sim,1),1]))...
            ./ repmat(precip_std, [size(precip_sim,1),1]);
        warning('Removing the seasonal cycle from occurrence data is not very intuitive. Are you sure this is what you want?');
    end
    
    % Shift the data so that it starts on the given start_DOY:
    precip_obs = ShiftXdays(precip_obs,1-start_DOY);
    precip_obs(end,:) = []; % This one has both the incomplete first year and the NaNs from the end of this year.
    precip_sim = ShiftXdays(precip_sim,1-start_DOY);
    precip_sim(end,:) = [];
    years(1) = [];

    % % % % % % % % % % % % % % %
    % 1yr
    if any(strcmpi('oneyear',periods)) || any(strcmpi('1yr',periods))

        precip_obs_tmp = precip_obs;
        precip_sim_tmp = precip_sim;

        % Now we need to make the sim data an integer multiple of the obs data
        % and insert NaNs at the correct places:
        n_years = size(precip_obs_tmp,1);
        n_sims = floor(size(precip_sim_tmp,1)/n_years);
        precip_sim_tmp( (n_years*n_sims+1):end, : ) = [];
        obs_nans = isnan(precip_obs_tmp);
        precip_sim_tmp( repmat(obs_nans,[n_sims,1]) ) = nan;
        
        precip_obs_1yr = 365*nanmean(precip_obs_tmp,2); % n_years x 1
        precip_sim_1yr = 365*reshape(nanmean(precip_sim_tmp,2), [n_years,n_sims]); % As an n_years x n_sims matrix
        
        % for occurrence:
        precip_obs_tmp(precip_obs_tmp>0) = 1; % Now occurrence!
        precip_sim_tmp(precip_sim_tmp>0) = 1; % Now occurrence!
        
        occ_obs_1yr = 365*nanmean(precip_obs_tmp,2); % n_years x 1
        occ_sim_1yr = 365*reshape(nanmean(precip_sim_tmp,2), [n_years,n_sims]); % As an n_years x n_sims matrix
        
		nino_djfma_1yr = nino_djfma;
		nino_djfma_1yr( 1:(length(nino_djfma_1yr)-length(occ_obs_1yr)) ) = [];
    end % if
    
    % % % % % % % % % % % % % % %
    % 2yr
    if any(strcmpi('twoyear',periods)) || any(strcmpi('2yr',periods))
        precip_obs_tmp = precip_obs;
        precip_sim_tmp = precip_sim;
        
        % Make the precip matrices have an even number of years:
        years_to_remove = mod(size(precip_obs,1),2);
        precip_obs_tmp(1:years_to_remove,:) = []; % Even!
        
        % Now we need to make the sim data an integer multiple of the obs data
        % and insert NaNs at the correct places:
        n_years = size(precip_obs_tmp,1);
        n_sims = floor(size(precip_sim_tmp,1)/n_years);
        precip_sim_tmp( (n_years*n_sims+1):end, : ) = [];
        obs_nans = isnan(precip_obs_tmp);
        precip_sim_tmp( repmat(obs_nans,[n_sims,1]) ) = nan;

        % Now make each row represent 2 years:
        precip_obs_tmp = reshape(precip_obs_tmp',[2*365, numel(precip_obs_tmp)/(2*365)])';
        precip_sim_tmp = reshape(precip_sim_tmp',[2*365, numel(precip_sim_tmp)/(2*365)])';
        
        precip_obs_2yr = 2*365*nanmean(precip_obs_tmp,2); % n_years x 1
        precip_sim_2yr = 2*365*reshape(nanmean(precip_sim_tmp,2),... 
            [length(precip_obs_2yr),size(precip_sim_tmp,1)/length(precip_obs_2yr)]); % As an n_years x n_sims matrix        

        % for occurrence:
        precip_obs_tmp(precip_obs_tmp>0) = 1; % Now occurrence!
        precip_sim_tmp(precip_sim_tmp>0) = 1; % Now occurrence!
        
        occ_obs_2yr = 2*365*nanmean(precip_obs_tmp,2); % n_years x 1
        occ_sim_2yr = 2*365*reshape(nanmean(precip_sim_tmp,2),... 
            [length(precip_obs_2yr),size(precip_sim_tmp,1)/length(precip_obs_2yr)]); % As an n_years x n_sims matrix        

		nino_djfma_2yr = nino_djfma;
		if mod(length(nino_djfma_2yr),2) == 1 % odd
			nino_djfma_2yr(1) = [];
		end
		% Combine neighboring years
		nino_djfma_2yr = reshape(nino_djfma_2yr,[2,length(nino_djfma_2yr)/2])';
		nino_djfma_2yr = sum(nino_djfma_2yr,2);
		nino_djfma_2yr( 1:(length(nino_djfma_2yr)-length(occ_obs_2yr)) ) = [];

    end % if for 2yr
    
    
    % % % % % % % % % % % % % % %
    % 3yr
    if any(strcmpi('threeyear',periods)) || any(strcmpi('3yr',periods))
        precip_obs_tmp = precip_obs;
        precip_sim_tmp = precip_sim;

        % Make the precip matrices have a number of years that is a multiple of 3:
        years_to_remove = mod(size(precip_obs,1),3);
        precip_obs_tmp(1:years_to_remove,:) = []; % Multiple of 3!
        
        % Now we need to make the sim data an integer multiple of the obs data
        % and insert NaNs at the correct places:
        n_years = size(precip_obs_tmp,1);
        n_sims = floor(size(precip_sim_tmp,1)/n_years);
        precip_sim_tmp( (n_years*n_sims+1):end, : ) = [];
        obs_nans = isnan(precip_obs_tmp);
        precip_sim_tmp( repmat(obs_nans,[n_sims,1]) ) = nan;

        % Now make each row represent 3 years:
        precip_obs_tmp = reshape(precip_obs_tmp',[3*365, numel(precip_obs_tmp)/(3*365)])';
        precip_sim_tmp = reshape(precip_sim_tmp',[3*365, numel(precip_sim_tmp)/(3*365)])';
     
        precip_obs_3yr = 3*365*nanmean(precip_obs_tmp,2); % n_years x 1
        precip_sim_3yr = 3*365*reshape(nanmean(precip_sim_tmp,2),... 
            [length(precip_obs_3yr),size(precip_sim_tmp,1)/length(precip_obs_3yr)]); % As an n_years x n_sims matrix        

        % for occurrence:
        precip_obs_tmp(precip_obs_tmp>0) = 1; % Now occurrence!
        precip_sim_tmp(precip_sim_tmp>0) = 1; % Now occurrence!

        occ_obs_3yr = 3*365*nanmean(precip_obs_tmp,2); % n_years x 1
        occ_sim_3yr = 3*365*reshape(nanmean(precip_sim_tmp,2),... 
            [length(precip_obs_3yr),size(precip_sim_tmp,1)/length(precip_obs_3yr)]); % As an n_years x n_sims matrix        

		nino_djfma_3yr = nino_djfma;
		if mod(length(nino_djfma_3yr),3) == 1
			nino_djfma_3yr(1) = [];
		elseif mod(length(nino_djfma_3yr),3) == 2
			nino_djfma_3yr(1:2) = [];
		end
		% Combine neighboring years
		nino_djfma_3yr = reshape(nino_djfma_3yr,[3,length(nino_djfma_3yr)/3])';
		nino_djfma_3yr = sum(nino_djfma_3yr,2);
		nino_djfma_3yr( 1:(length(nino_djfma_3yr)-length(occ_obs_3yr)) ) = [];


    end % if for 3yr

    % Now save the new precip time series:
    filename = sprintf('precip_%s_accum_DOY%d.mat',id,start_DOY);
    save(filename, 'precip_obs_1yr','precip_sim_1yr', 'precip_obs_2yr',...
        'precip_sim_2yr','precip_obs_3yr','precip_sim_3yr',...
        'occ_obs_1yr','occ_sim_1yr', 'occ_obs_2yr','occ_sim_2yr',...
        'occ_obs_3yr','occ_sim_3yr','id','start_DOY',...
		'nino_djfma_1yr','nino_djfma_2yr','nino_djfma_3yr');
    
end




end