function std_norm_pdfs_obs = sim_dist_to_obs_std_norm(LL_obs, LL_sim, uniform_pdfs_sim)

X = [LL_sim(:), uniform_pdfs_sim(:)];
sorted = sortrows(X, 1); % Sort by LL_sim values (sorting by column 2 should give the same ans


% Now, we just interpolate the points that are within the LL_sim
% distribution:
uniform_pdfs_obs = interp1(sorted(:,1),sorted(:,2),LL_obs);

% But we have to extrapolate for the LL_obs values outside that distribution:
extrapolated_points = isnan(uniform_pdfs_obs);

% here's the plan: Assume we have a point with a larger LL_obs than any in
% LL_sim -- we fit a normal distribution to LL_sim, find the cdf value of
% the largest sim value, find the cdf value of the obs point in that normal
% distribution (it will be higher), calculate the ratio of the cdfs, and
% then assign the uniform pdf as the same ratio of the max of the uniform
% sim pdfs:
if(sum(extrapolated_points) > 0)
    fprintf('We are extrapolating obs LL values! Too bad. =(\n');

    [mu,sigma] = normfit(LL_sim);
    max_sim = max(LL_sim(:));
    min_sim = min(LL_sim(:));
    
    high_vals = LL_obs > max_sim;
    low_vals = LL_obs < min_sim;
    
    sim_max_cdf = normcdf(max_sim,mu,sigma);
    sim_min_cdf = normcdf(min_sim,mu,sigma);
    
    high_ratios = normcdf( LL_obs(high_vals), mu, sigma ) / sim_max_cdf;
    low_ratios = normcdf( LL_obs(low_vals), mu, sigma ) / sim_min_cdf;

    uniform_pdfs_obs(high_vals) = max(uniform_pdfs_sim(:)) * high_ratios; % should all be larger than max(sim_uniform), but less than 1
    uniform_pdfs_obs(low_vals) = min(uniform_pdfs_sim(:)) * low_ratios; % should all be smaller than min(sim_uniform), but larger than 0    
end


std_norm_pdfs_obs = norminv(uniform_pdfs_obs,0,1); 


end % function