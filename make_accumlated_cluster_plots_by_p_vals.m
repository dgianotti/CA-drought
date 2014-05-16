function make_accumlated_cluster_plots_by_p_vals(DOY)

load('CA_ids.mat');

marker_size = 40;
clim = [-5,5];

%load('cmap_better_jet');
%cmap = flipud(cmap_better_jet);

% Or, better yet,
addpath('C:\Users\gianotti\Documents\GitHub\CA-drought\cbrewer\cbrewer\');
cmap = cbrewer('div','RdYlBu',5,'cubic');
n = floor(size(cmap,1)/2);
cmap = [cmap(1:n,:); cmap(n+1,:); cmap(n+1,:); cmap(n+1,:); cmap( (n+2):end,:)];

for i = 1:length(good_CA_IDs)
    fprintf('Clustering station %i of %i...\n',i,length(good_CA_IDs));
    stn_id = good_CA_IDs{i};
    % Load the data:
    filename = sprintf('LL_%s_accum_DOY%i.mat',stn_id,DOY);    
    load(filename);

    filename = sprintf('precip_%s_accum_DOY%i.mat',stn_id,DOY);
    load(filename);
    
    % % %% % %% % % %% 
    % 1-yearly data:
    
    % Cluster the LL_obs data:
    [clusters, mu_hat, sigma_hat, weights] = GMM_cluster(LL_obs_1yr_stdnorm, LL_sim_1yr_stdnorm);

    N = length(LL_obs_1yr_stdnorm);
    years = (2014-N):2013;
    if (DOY > 1) % using something later than Jan1 start date
        years = years+1; % so that it ends with 2014
    end
    
    figure;
    subplot(3,6,1:4);
    plot(years, -LL_obs_1yr_stdnorm,'-k');
    hold on;
    scatter(years, -LL_obs_1yr_stdnorm,marker_size,-LL_obs_1yr_stdnorm,'filled','MarkerEdgeColor','k');
    colormap(cmap);
    caxis(clim);
    
    ylm = get(gca,'Ylim');
    title(sprintf('Stn #%s: 1-yearly, Start DOY = %i',stn_id,DOY));
    ylabel('-Log-Likelihood');
    
    subplot(3,6,5);
    y = linspace(ylm(1),ylm(2),200);
    plot(weights(1)*normpdf(y,mu_hat(1),sigma_hat(1)),-y,'-b');
    hold on;
    plot(weights(2)*normpdf(y,mu_hat(2),sigma_hat(2)),-y,'-k');
    plot(weights(3)*normpdf(y,mu_hat(3),sigma_hat(3)),-y,'-r');
    ylim(ylm);
    
    % QQ plot:
    subplot(3,6,6);
    qqplot(-LL_sim_1yr_stdnorm(:),-LL_obs_1yr_stdnorm(:));
    xlabel('Sim Quantile');
    ylabel('Obs Quantile');
    ylim(ylm);

    % % % % % % % % % 
    % Now the precip data:
    [clusters, mu_hat, sigma_hat, weights] = GMM_cluster(precip_obs_1yr_stdnorm, precip_sim_1yr_stdnorm);
    subplot(3,6,7:10);
    plot(years, precip_obs_1yr_stdnorm,'-k');
    hold on;
    scatter(years, precip_obs_1yr_stdnorm,marker_size,precip_obs_1yr_stdnorm,'filled','MarkerEdgeColor','k');
    colormap(cmap);
    caxis(clim);
    
    ylm = get(gca,'Ylim');
    ylabel('Total Precipitation');
    
    subplot(3,6,11);
    y = linspace(ylm(1),ylm(2),200);
    plot(weights(1)*normpdf(y,mu_hat(1),sigma_hat(1)),y,'-r');
    hold on;
    plot(weights(2)*normpdf(y,mu_hat(2),sigma_hat(2)),y,'-k');
    plot(weights(3)*normpdf(y,mu_hat(3),sigma_hat(3)),y,'-b');
    ylim(ylm);

    % QQ plot:
    subplot(3,6,12);
    qqplot(precip_sim_1yr_stdnorm(:),precip_obs_1yr_stdnorm(:));
    xlabel('Sim Quantile');
    ylabel('Obs Quantile');
    ylim(ylm);
    
    % % % % % % % % % 
    % Now the occ data:
    [clusters, mu_hat, sigma_hat, weights] = GMM_cluster(occ_obs_1yr_stdnorm, occ_sim_1yr_stdnorm);
    subplot(3,6,13:16);
    plot(years, occ_obs_1yr_stdnorm,'-k');
    hold on;
    scatter(years, occ_obs_1yr_stdnorm,marker_size,occ_obs_1yr_stdnorm,'filled','MarkerEdgeColor','k');
    colormap(cmap);
    caxis(clim);

    ylm = get(gca,'Ylim');
    ylabel({'Precipitation','Occurrence'});
    
    subplot(3,6,17);
    y = linspace(ylm(1),ylm(2),200);
    plot(weights(1)*normpdf(y,mu_hat(1),sigma_hat(1)),y,'-r');
    hold on;
    plot(weights(2)*normpdf(y,mu_hat(2),sigma_hat(2)),y,'-k');
    plot(weights(3)*normpdf(y,mu_hat(3),sigma_hat(3)),y,'-b');
    ylim(ylm);
    
    % QQ plot:
    subplot(3,6,18);
    qqplot(occ_sim_1yr_stdnorm(:),occ_obs_1yr_stdnorm(:));
    xlabel('Sim Quantile');
    ylabel('Obs Quantile');
    ylim(ylm);
    
    outfilename = sprintf('plots/Pval_TS_%s_1yr_DOY%i.png',stn_id,DOY);
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
    subplot(3,6,1:4);
    plot(years, -LL_obs_2yr_stdnorm,'-k');
    hold on;
    scatter(years, -LL_obs_2yr_stdnorm,marker_size,-LL_obs_2yr_stdnorm,'filled','MarkerEdgeColor','k');
    colormap(cmap);
    caxis(clim);

    ylm = get(gca,'Ylim');
    title(sprintf('Stn #%s: 2-yearly, Start DOY = %i',stn_id,DOY));
    ylabel('-Log-Likelihood');
    
    subplot(3,6,5);
    y = linspace(ylm(1),ylm(2),200);
    plot(weights(1)*normpdf(y,mu_hat(1),sigma_hat(1)),-y,'-b');
    hold on;
    plot(weights(2)*normpdf(y,mu_hat(2),sigma_hat(2)),-y,'-k');
    plot(weights(3)*normpdf(y,mu_hat(3),sigma_hat(3)),-y,'-r');
    ylim(ylm);
    
    % QQ plot:
    subplot(3,6,6);
    qqplot(-LL_sim_2yr_stdnorm(:),-LL_obs_2yr_stdnorm(:));
    xlabel('Sim Quantile');
    ylabel('Obs Quantile');
    ylim(ylm);

    % % % % % % % % % 
    % Now the precip data:
    [clusters, mu_hat, sigma_hat, weights] = GMM_cluster(precip_obs_2yr_stdnorm, precip_sim_2yr_stdnorm);
    subplot(3,6,7:10);
    plot(years, precip_obs_2yr_stdnorm,'-k');
    hold on;
    scatter(years, precip_obs_2yr_stdnorm,marker_size,precip_obs_2yr_stdnorm,'filled','MarkerEdgeColor','k');
    colormap(cmap);
    caxis(clim);

    ylm = get(gca,'Ylim');
    ylabel('Total Precipitation');
    
    subplot(3,6,11);
    y = linspace(ylm(1),ylm(2),200);
    plot(weights(1)*normpdf(y,mu_hat(1),sigma_hat(1)),y,'-r');
    hold on;
    plot(weights(2)*normpdf(y,mu_hat(2),sigma_hat(2)),y,'-k');
    plot(weights(3)*normpdf(y,mu_hat(3),sigma_hat(3)),y,'-b');
    ylim(ylm);

    % QQ plot:
    subplot(3,6,12);
    qqplot(precip_sim_2yr_stdnorm(:),precip_obs_2yr_stdnorm(:));
    xlabel('Sim Quantile');
    ylabel('Obs Quantile');
    ylim(ylm);
    
    % % % % % % % % % 
    % Now the occ data:
    [clusters, mu_hat, sigma_hat, weights] = GMM_cluster(occ_obs_2yr_stdnorm, occ_sim_2yr_stdnorm);
    subplot(3,6,13:16);
    plot(years, occ_obs_2yr_stdnorm,'-k');
    hold on;
    scatter(years, occ_obs_2yr_stdnorm,marker_size,occ_obs_2yr_stdnorm,'filled','MarkerEdgeColor','k');
    colormap(cmap);
    caxis(clim);

    ylm = get(gca,'Ylim');
    ylabel({'Precipitation','Occurrence'});
    
    subplot(3,6,17);
    y = linspace(ylm(1),ylm(2),200);
    plot(weights(1)*normpdf(y,mu_hat(1),sigma_hat(1)),y,'-r');
    hold on;
    plot(weights(2)*normpdf(y,mu_hat(2),sigma_hat(2)),y,'-k');
    plot(weights(3)*normpdf(y,mu_hat(3),sigma_hat(3)),y,'-b');
    ylim(ylm);
    
    
    % QQ plot:
    subplot(3,6,18);
    qqplot(occ_sim_2yr_stdnorm(:),occ_obs_2yr_stdnorm(:));
    xlabel('Sim Quantile');
    ylabel('Obs Quantile');
    ylim(ylm);
    
    outfilename = sprintf('plots/Pval_TS_%s_2yr_DOY%i.png',stn_id,DOY);
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
    subplot(3,6,1:4);
    plot(years, -LL_obs_3yr_stdnorm,'-k');
    hold on;
    scatter(years, -LL_obs_3yr_stdnorm,marker_size,-LL_obs_3yr_stdnorm,'filled','MarkerEdgeColor','k');
    colormap(cmap);
    caxis(clim);

    ylm = get(gca,'Ylim');
    title(sprintf('Stn #%s: 3-yearly, Start DOY = %i',stn_id,DOY));
    ylabel('-Log-Likelihood');
   
    subplot(3,6,5);
    y = linspace(ylm(1),ylm(2),200);
    plot(weights(1)*normpdf(y,mu_hat(1),sigma_hat(1)),-y,'-b');
    hold on;
    plot(weights(2)*normpdf(y,mu_hat(2),sigma_hat(2)),-y,'-k');
    plot(weights(3)*normpdf(y,mu_hat(3),sigma_hat(3)),-y,'-r');
    ylim(ylm);

    % QQ plot:
    subplot(3,6,6);
    qqplot(-LL_sim_3yr_stdnorm(:),-LL_obs_3yr_stdnorm(:));
    xlabel('Sim Quantile');
    ylabel('Obs Quantile');
    ylim(ylm);

    % % % % % % % % % 
    % Now the precip data:
    [clusters, mu_hat, sigma_hat, weights] = GMM_cluster(precip_obs_3yr_stdnorm, precip_sim_3yr_stdnorm);
    subplot(3,6,7:10);
    plot(years, precip_obs_3yr_stdnorm,'-k');
    hold on;
    scatter(years, precip_obs_3yr_stdnorm,marker_size,precip_obs_3yr_stdnorm,'filled','MarkerEdgeColor','k');
    colormap(cmap);
    caxis(clim);

    ylm = get(gca,'Ylim');
    ylabel('Total Precipitation');
    
    subplot(3,6,11);
    y = linspace(ylm(1),ylm(2),200);
    plot(weights(1)*normpdf(y,mu_hat(1),sigma_hat(1)),y,'-r');
    hold on;
    plot(weights(2)*normpdf(y,mu_hat(2),sigma_hat(2)),y,'-k');
    plot(weights(3)*normpdf(y,mu_hat(3),sigma_hat(3)),y,'-b');
    ylim(ylm);

    % QQ plot:
    subplot(3,6,12);
    qqplot(precip_sim_3yr_stdnorm(:),precip_obs_3yr_stdnorm(:));
    xlabel('Sim Quantile');
    ylabel('Obs Quantile');
    ylim(ylm);
    
    % % % % % % % % % 
    % Now the occ data:
    [clusters, mu_hat, sigma_hat, weights] = GMM_cluster(occ_obs_3yr_stdnorm, occ_sim_3yr_stdnorm);
    subplot(3,6,13:16);
    plot(years, occ_obs_3yr_stdnorm,'-k');
    hold on;
    scatter(years, occ_obs_3yr_stdnorm,marker_size,occ_obs_3yr_stdnorm,'filled','MarkerEdgeColor','k');
    colormap(cmap);
    caxis(clim);

    ylm = get(gca,'Ylim');
    ylabel({'Precipitation','Occurrence'});
    
    subplot(3,6,17);
    y = linspace(ylm(1),ylm(2),200);
    plot(weights(1)*normpdf(y,mu_hat(1),sigma_hat(1)),y,'-r');
    hold on;
    plot(weights(2)*normpdf(y,mu_hat(2),sigma_hat(2)),y,'-k');
    plot(weights(3)*normpdf(y,mu_hat(3),sigma_hat(3)),y,'-b');
    ylim(ylm);
    
    % QQ plot:
    subplot(3,6,18);
    qqplot(occ_sim_3yr_stdnorm(:),occ_obs_3yr_stdnorm(:));
    xlabel('Sim Quantile');
    ylabel('Obs Quantile');
    ylim(ylm);
    
    outfilename = sprintf('plots/Pval_TS_%s_3yr_DOY%i.png',stn_id,DOY);
    print(gcf,'-dpng',outfilename);
    
end % for loop over stations
    
end % function