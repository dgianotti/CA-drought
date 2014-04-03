function std_norm_pdfs_obs = sim_dist_to_obs_std_norm(LL_obs, LL_sim, uniform_pdfs_sim)

X = [LL_sim(:), uniform_pdfs_sim(:)];
sorted = sortrows(X, 1); % Sort by LL_sim values (sorting by column 2 should give the same ans


% Now, we just interpolate the points that are within the LL_sim
% distribution:
uniform_pdfs_obs = interp1(sorted(:,1),sorted(:,2),LL_obs);

% But we have to extrapolate for the LL_obs values outside that distribution:
extrapolated_points = isnan(uniform_pdfs_obs);


std_norm_pdfs_obs = norminv(uniform_pdfs_obs,0,1); 


end % function