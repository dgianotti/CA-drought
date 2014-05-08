function cdf_transform_accum_precip_data()

% Get the file names:
accum_files = dir('precip*accum_DOY*.mat');

% Number of files
num_files = numel(accum_files);

for i = 1:num_files
    fprintf('Transforming file %i of %i...\n',i,num_files);
    % Load the data:
    load(accum_files(i).name);
    
    precip_obs_1yr_stdnorm = norminv( ksdensity(precip_sim_1yr(:),precip_obs_1yr,'function','cdf'), 0, 1);
    precip_sim_1yr_stdnorm = (precip_sim_1yr-mean(precip_sim_1yr(:))) ./ std(precip_sim_1yr(:));
 
    precip_obs_2yr_stdnorm = norminv( ksdensity(precip_sim_2yr(:),precip_obs_2yr,'function','cdf'), 0, 1);
    precip_sim_2yr_stdnorm = (precip_sim_2yr-mean(precip_sim_2yr(:))) ./ std(precip_sim_2yr(:));
 
    precip_obs_3yr_stdnorm = norminv( ksdensity(precip_sim_3yr(:),precip_obs_3yr,'function','cdf'), 0, 1);
    precip_sim_3yr_stdnorm = (precip_sim_3yr-mean(precip_sim_3yr(:))) ./ std(precip_sim_3yr(:)) ;

    occ_obs_1yr_stdnorm = norminv( ksdensity(occ_sim_1yr(:),occ_obs_1yr,'function','cdf'), 0, 1);
    occ_sim_1yr_stdnorm = (occ_sim_1yr-mean(occ_sim_1yr(:))) ./ std(occ_sim_1yr(:));
 
    occ_obs_2yr_stdnorm = norminv( ksdensity(occ_sim_2yr(:),occ_obs_2yr,'function','cdf'), 0, 1);
    occ_sim_2yr_stdnorm = (occ_sim_2yr-mean(occ_sim_2yr(:))) ./ std(occ_sim_2yr(:));
 
    occ_obs_3yr_stdnorm = norminv( ksdensity(occ_sim_3yr(:),occ_obs_3yr,'function','cdf'), 0, 1);
    occ_sim_3yr_stdnorm = (occ_sim_3yr-mean(occ_sim_3yr(:))) ./ std(occ_sim_3yr(:)) ;
    
    
    save(accum_files(i).name, 'precip_obs_1yr_stdnorm','precip_sim_1yr_stdnorm',...
        'precip_obs_2yr_stdnorm','precip_sim_2yr_stdnorm',...
        'precip_obs_3yr_stdnorm','precip_sim_3yr_stdnorm',...
        'occ_obs_1yr_stdnorm','occ_sim_1yr_stdnorm',...
        'occ_obs_2yr_stdnorm','occ_sim_2yr_stdnorm',...
        'occ_obs_3yr_stdnorm','occ_sim_3yr_stdnorm','-append');
end


end