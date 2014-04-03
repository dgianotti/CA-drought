% Try the same thing as in CompileStnData, but with synthetic data:

clear;
clc;

if strcmpi(getenv('OS'),'Windows_NT')
    addpath('C:\Users\gianotti\Documents\IntensityLib\')
else % linux
    addpath('IntensityLib');
end

addpath('..');

today_doy = 1; % Jan1

id = 'synthetic';

% Make the daily precip data:
% obs_params_1 = [0.3, .5, 3.5]; % p_occ, gamma_1, gamma_2
% obs_params_2 = [0.1, 2,2];

obs_params_1 = [0.3, 5, 2]; % p_occ, gamma_1, gamma_2
obs_params_2 = [0.1, 2, 15];

    

obs1 = gamrnd(obs_params_1(2),obs_params_1(3),[60,365]) .* (rand(60,365) < obs_params_1(1));
obs2 = gamrnd(obs_params_2(2),obs_params_2(3),[40,365]) .* (rand(40,365) < obs_params_2(1));
obs_precip = [obs1;obs2];

% Get synthetic model params:
p_occ = sum(obs_precip(:)>0)/numel(obs_precip)
obs_vec = reshape(obs_precip',numel(obs_precip),1);
gam_params = gamfit(obs_vec(obs_vec>0))

sim_precip = gamrnd(gam_params(1),gam_params(2),[100*1000,365]) .* (rand(100*1000,365) < p_occ);
sim_vec = reshape(sim_precip',numel(sim_precip),1);

% Calculate likelihoods:
%% LL_obs
cdfs = zeros(size(obs_vec));
cdfs(obs_vec == 0) = 1 - p_occ;
cdfs(obs_vec > 0) = 1 - p_occ + p_occ* gamcdf(obs_vec(obs_vec > 0),gam_params(1),gam_params(2));
LL_obs_std_norm = log(normpdf(norminv(cdfs,0,1),0,1)); % norminv(cdfs) is a precip equivalent, normlike(norminv(cdfs)) is -LL
LL_obs_std_norm = reshape(LL_obs_std_norm,[365,numel(LL_obs_std_norm)/365])';

%% LL_sim
cdfs = zeros(size(sim_vec));
cdfs(sim_vec == 0) = 1 - p_occ;
cdfs(sim_vec > 0) = 1 - p_occ + p_occ*gamcdf(sim_vec(sim_vec > 0),gam_params(1),gam_params(2));
LL_sim_std_norm = log(normpdf(norminv(cdfs,0,1),0,1));
LL_sim_std_norm = reshape(LL_sim_std_norm,[365,numel(LL_sim_std_norm)/365])';
    

%%

%load(['LL_',id,'.mat']); 

% Now we have LL_obs, LL_sim, LL_obs_std_norm, LL_sim_std_norm, and years

% If you want to use the std_norm versions, uncomment the following:
LL_obs = LL_obs_std_norm;
LL_sim = LL_sim_std_norm;
clear LL_sim_std_norm LL_obs_std_norm;

% Remove the seasonal cycle as best you can:
LL_mean = mean(LL_sim,1);
LL_std = std(LL_sim,0,1);
LL_obs = (LL_obs - repmat(LL_mean, [size(LL_obs,1),1]))... 
    ./ repmat(LL_std, [size(LL_obs,1),1]);
LL_sim = (LL_sim - repmat(LL_mean, [size(LL_sim,1),1]))... 
    ./ repmat(LL_std, [size(LL_sim,1),1]);

% Okay, now they're a little more Gaussian...

% Calculate the annual LL, with the year begining on today's DOY:
LL_obs_shifted = ShiftXdays(LL_obs,1-today_doy);
LL_sim_shifted = ShiftXdays(LL_sim,1-today_doy);

% Now we have to turn the sim data into an integer multiple of the obs
% data, and make sure that it has nans in the same place:
n_years = size(LL_obs_shifted,1);
n_sims = floor(size(LL_sim,1)/n_years);
LL_sim_shifted( (n_years*n_sims+1):end, : ) = [];
obs_nans = isnan(LL_obs_shifted);
LL_sim_shifted( repmat(obs_nans,[n_sims,1]) ) = nan;

LL_obs_annual = nanmean(LL_obs_shifted,2); % n_years x 1
LL_sim_annual = reshape(nanmean(LL_sim_shifted,2), [n_years,n_sims]); % As an n_years x n_sims matrix

% Convert the sim to a normal distribution using order statistics:
[std_norm_pdfs_sim, uniform_pdfs_sim, ~] = empirical_2_normal_via_order_stats(LL_sim_annual);

% Get the std_normal values for the obs using the sim transform as a
% look-up table:
std_norm_pdfs_obs = sim_dist_to_obs_std_norm(LL_obs_annual, LL_sim_annual, ...
    uniform_pdfs_sim);

years = 1901:2000;
plot(years,quantile(std_norm_pdfs_sim,[.025,.5,.975],2),'-r');
hold on;
plot(years,std_norm_pdfs_obs,'-b');
xlim([10*floor(min(years/10)),2020]);
title(id);
set(gca,'XTick',1900:50:2000);

print(gcf,'-dpdf','synthetic_AnnualNorm.pdf');

figure;
plot(years,quantile(LL_sim_annual,[.025,.5,.975],2),'-r');
hold on;
plot(years,LL_obs_annual,'-b');
xlim([10*floor(min(years/10)),2020]);
title(id);
set(gca,'XTick',1900:50:2000);

print(gcf,'-dpdf','synthetic_DailyNorm_NoAnnualNorm.pdf');
