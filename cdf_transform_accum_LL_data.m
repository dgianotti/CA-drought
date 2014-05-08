function cdf_transform_accum_LL_data()

% Get the file names:
accum_files = dir('LL*accum_DOY*.mat');

% Number of files
num_files = numel(accum_files);

for i = 1:num_files
    fprintf('Transforming file %i of %i...\n',i,num_files);
    % Load the data:
    load(accum_files(i).name);
    
    LL_obs_1yr_stdnorm = norminv( ksdensity(LL_sim_1yr(:),LL_obs_1yr,'function','cdf'), 0, 1);
    LL_sim_1yr_stdnorm = (LL_sim_1yr-mean(LL_sim_1yr(:))) ./ std(LL_sim_1yr(:)) ;
 
    LL_obs_2yr_stdnorm = norminv( ksdensity(LL_sim_2yr(:),LL_obs_2yr,'function','cdf'), 0, 1);
    LL_sim_2yr_stdnorm = (LL_sim_1yr-mean(LL_sim_2yr(:))) ./ std(LL_sim_2yr(:)) ;
 
    LL_obs_3yr_stdnorm = norminv( ksdensity(LL_sim_3yr(:),LL_obs_3yr,'function','cdf'), 0, 1);
    LL_sim_3yr_stdnorm = (LL_sim_1yr-mean(LL_sim_3yr(:))) ./ std(LL_sim_3yr(:)) ;

    save(accum_files(i).name, 'LL_obs_1yr_stdnorm','LL_sim_1yr_stdnorm',...
        'LL_obs_2yr_stdnorm','LL_sim_2yr_stdnorm',...
        'LL_obs_3yr_stdnorm','LL_sim_3yr_stdnorm', '-append');
end


end