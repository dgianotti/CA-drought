function make_accumlated_cluster_plots_LL(DOY)

load('CA_ids.mat');

for i = 1:length(good_CA_IDs)
    fprintf('Clustering station %i of %i...\n',i,length(good_CA_IDs));
    stn_id = good_CA_IDs{i};
    % Load the data:
    filename = sprintf('LL_%s_accum_DOY%i.mat',stn_id,DOY);
    
    load(filename);
    
    % % %% % %% % % %% 
    % 1-yearly data:
    
    % Cluster the obs data:
    [clusters, mu_hat, sigma_hat, weights] = GMM_clusterEM(LL_obs_1yr_stdnorm, LL_sim_1yr_stdnorm);

    N = length(LL_obs_1yr_stdnorm);
    years = (2014-N):2013;
    if (DOY > 1) % using something later than Jan1 start date
        years = years+1; % so that it ends with 2014
    end
    
    figure;
    subplot(1,5,1:4);
    plot(years(clusters == 1), LL_obs_1yr_stdnorm(clusters == 1),'ob',...
        'LineStyle','none');
    hold on;
    plot(years(clusters == 2), LL_obs_1yr_stdnorm(clusters == 2),'ok','LineStyle','none');
    plot(years(clusters == 3), LL_obs_1yr_stdnorm(clusters == 3),'or','LineStyle','none');
    ylm = get(gca,'Ylim');
    title(sprintf('Stn #%s: 1-yearly, Start DOY = %i',stn_id,DOY));
    
    subplot(1,5,5);
    y = linspace(ylm(1),ylm(2),200);
    plot(weights(1)*normpdf(y,mu_hat(1),sigma_hat(1)),y,'-b');
    hold on;
    plot(weights(2)*normpdf(y,mu_hat(2),sigma_hat(2)),y,'-k');
    plot(weights(3)*normpdf(y,mu_hat(3),sigma_hat(3)),y,'-r');
    ylim(ylm);
    % set(gca,'Xscale','log');
    
    outfilename = sprintf('plots/LL_TS_%s_1yr_DOY%i.png',stn_id,DOY);
    print(gcf,'-dpng',outfilename);
    
    
    % % %% % %% % % %% 
    % 2-yearly data:
    
    % Cluster the obs data:
    [clusters, mu_hat, sigma_hat, weights] = GMM_cluster(LL_obs_2yr_stdnorm, LL_sim_2yr_stdnorm);

    N = length(LL_obs_2yr_stdnorm);
    years = (2014.5-2*N):2:2012.5;
    if (DOY > 1) % using something later than Jan1 start date
        years = years+1; % so that it ends with 2014
    end
    
    figure;
    subplot(1,5,1:4);
    plot(years(clusters == 1), LL_obs_2yr_stdnorm(clusters == 1),'ob',...
        'LineStyle','none');
    hold on;
    plot(years(clusters == 2), LL_obs_2yr_stdnorm(clusters == 2),'ok','LineStyle','none');
    plot(years(clusters == 3), LL_obs_2yr_stdnorm(clusters == 3),'or','LineStyle','none');
    ylm = get(gca,'Ylim');
    title(sprintf('Stn #%s: 2-yearly, Start DOY = %i',stn_id,DOY));
    
    subplot(1,5,5);
    y = linspace(ylm(1),ylm(2),200);
    plot(weights(1)*normpdf(y,mu_hat(1),sigma_hat(1)),y,'-b');
    hold on;
    plot(weights(2)*normpdf(y,mu_hat(2),sigma_hat(2)),y,'-k');
    plot(weights(3)*normpdf(y,mu_hat(3),sigma_hat(3)),y,'-r');
    ylim(ylm);
    
    outfilename = sprintf('plots/LL_TS_%s_2yr_DOY%i.png',stn_id,DOY);
    print(gcf,'-dpng',outfilename);
    
    
    % % %% % %% % % %%
    % 3-yearly data:
    
    % Cluster the obs data:
    [clusters, mu_hat, sigma_hat, weights] = GMM_cluster(LL_obs_3yr_stdnorm, LL_sim_3yr_stdnorm);

    N = length(LL_obs_3yr_stdnorm);
    years = (2015-3*N):3:2012;
    if (DOY > 1) % using something later than Jan1 start date
        years = years+1; % so that it ends with 2014
    end
    
    figure;
    subplot(1,5,1:4);
    plot(years(clusters == 1), LL_obs_3yr_stdnorm(clusters == 1),'ob',...
        'LineStyle','none');
    hold on;
    plot(years(clusters == 2), LL_obs_3yr_stdnorm(clusters == 2),'ok','LineStyle','none');
    plot(years(clusters == 3), LL_obs_3yr_stdnorm(clusters == 3),'or','LineStyle','none');
    ylm = get(gca,'Ylim');
    title(sprintf('Stn #%s: 3-yearly, Start DOY = %i',stn_id,DOY));
    
    subplot(1,5,5);
    y = linspace(ylm(1),ylm(2),200);
    plot(weights(1)*normpdf(y,mu_hat(1),sigma_hat(1)),y,'-b');
    hold on;
    plot(weights(2)*normpdf(y,mu_hat(2),sigma_hat(2)),y,'-k');
    plot(weights(3)*normpdf(y,mu_hat(3),sigma_hat(3)),y,'-r');
    ylim(ylm);
    
    outfilename = sprintf('plots/LL_TS_%s_3yr_DOY%i.png',stn_id,DOY);
    print(gcf,'-dpng',outfilename);
    
end % for loop over stations
    
end % function