function make_accumulated_cluster_plots_wet_season(DOY)

load('CA_ids.mat');

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
    subplot(7,5,1:4);
    plot(years, -LL_obs_1yr_stdnorm,'-k');
    hold on;
    plot(years(clusters == 1), -LL_obs_1yr_stdnorm(clusters == 1),'ob',...
        'LineStyle','none','MarkerFaceColor','b');
    plot(years(clusters == 2), -LL_obs_1yr_stdnorm(clusters == 2),'ok',...
        'LineStyle','none','MarkerFaceColor','k');
    plot(years(clusters == 3), -LL_obs_1yr_stdnorm(clusters == 3),'or',...
        'LineStyle','none','MarkerFaceColor','r');
    ylm = get(gca,'Ylim');
    title(sprintf('Stn #%s: 1-yearly, Start DOY = %i',stn_id,DOY));
    ylabel('-Log-Likelihood');
    
    subplot(7,5,5);
    y = linspace(ylm(1),ylm(2),200);
    plot(weights(1)*normpdf(y,mu_hat(1),sigma_hat(1)),-y,'-b');
    hold on;
    plot(weights(2)*normpdf(y,mu_hat(2),sigma_hat(2)),-y,'-k');
    plot(weights(3)*normpdf(y,mu_hat(3),sigma_hat(3)),-y,'-r');
    ylim(ylm);

    % % % % % % % % % 
    % Now the precip data:
    [clusters, mu_hat, sigma_hat, weights] = GMM_cluster(precip_obs_1yr_stdnorm, precip_sim_1yr_stdnorm);
    subplot(7,5,6:9);
    plot(years, precip_obs_1yr_stdnorm,'-k');
    hold on;
    plot(years(clusters == 1), precip_obs_1yr_stdnorm(clusters == 1),'or',...
        'LineStyle','none','MarkerFaceColor','r');
    plot(years(clusters == 2), precip_obs_1yr_stdnorm(clusters == 2),'ok',...
        'LineStyle','none','MarkerFaceColor','k');
    plot(years(clusters == 3), precip_obs_1yr_stdnorm(clusters == 3),'ob',...
        'LineStyle','none','MarkerFaceColor','b');
    ylm = get(gca,'Ylim');
    ylabel('Total Precipitation');
    
    subplot(7,5,10);
    y = linspace(ylm(1),ylm(2),200);
    plot(weights(1)*normpdf(y,mu_hat(1),sigma_hat(1)),y,'-r');
    hold on;
    plot(weights(2)*normpdf(y,mu_hat(2),sigma_hat(2)),y,'-k');
    plot(weights(3)*normpdf(y,mu_hat(3),sigma_hat(3)),y,'-b');
    ylim(ylm);

    % % % % % % % % % 
    % Now the occ data:
    [clusters, mu_hat, sigma_hat, weights] = GMM_cluster(occ_obs_1yr_stdnorm, occ_sim_1yr_stdnorm);
    subplot(7,5,11:14);
    plot(years, occ_obs_1yr_stdnorm,'-k');
    hold on;
    plot(years(clusters == 1), occ_obs_1yr_stdnorm(clusters == 1),'or',...
        'LineStyle','none','MarkerFaceColor','r');
    plot(years(clusters == 2), occ_obs_1yr_stdnorm(clusters == 2),'ok',...
        'LineStyle','none','MarkerFaceColor','k');
    plot(years(clusters == 3), occ_obs_1yr_stdnorm(clusters == 3),'ob',...
        'LineStyle','none','MarkerFaceColor','b');
    ylm = get(gca,'Ylim');
    ylabel({'Precipitation','Occurrence'});
    
    subplot(7,5,15);
    y = linspace(ylm(1),ylm(2),200);
    plot(weights(1)*normpdf(y,mu_hat(1),sigma_hat(1)),y,'-r');
    hold on;
    plot(weights(2)*normpdf(y,mu_hat(2),sigma_hat(2)),y,'-k');
    plot(weights(3)*normpdf(y,mu_hat(3),sigma_hat(3)),y,'-b');
    ylim(ylm);
    
    
    %% Now for the winter only precip:
    filename = sprintf('LL_%s_accum_DOY%i_wet_season_only.mat',stn_id,DOY);
    load(filename);

    filename = sprintf('precip_%s_accum_DOY%i_wet_season_only.mat',stn_id,DOY);
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

    subplot(7,5,20+(1:4));
    plot(years, -LL_obs_1yr_stdnorm,'-k');
    hold on;
    plot(years(clusters == 1), -LL_obs_1yr_stdnorm(clusters == 1),'ob',...
        'LineStyle','none','MarkerFaceColor','b');
    plot(years(clusters == 2), -LL_obs_1yr_stdnorm(clusters == 2),'ok',...
        'LineStyle','none','MarkerFaceColor','k');
    plot(years(clusters == 3), -LL_obs_1yr_stdnorm(clusters == 3),'or',...
        'LineStyle','none','MarkerFaceColor','r');
    ylm = get(gca,'Ylim');
    title(sprintf('Stn #%s: 1-yearly, Start DOY = %i -- WET SEASON ONLY',stn_id,DOY));
    ylabel('-Log-Likelihood');
    
    subplot(7,5,20+5);
    y = linspace(ylm(1),ylm(2),200);
    plot(weights(1)*normpdf(y,mu_hat(1),sigma_hat(1)),-y,'-b');
    hold on;
    plot(weights(2)*normpdf(y,mu_hat(2),sigma_hat(2)),-y,'-k');
    plot(weights(3)*normpdf(y,mu_hat(3),sigma_hat(3)),-y,'-r');
    ylim(ylm);

    % % % % % % % % % 
    % Now the precip data:
    [clusters, mu_hat, sigma_hat, weights] = GMM_cluster(precip_obs_1yr_stdnorm, precip_sim_1yr_stdnorm);
    subplot(7,5,20+(6:9));
    plot(years, precip_obs_1yr_stdnorm,'-k');
    hold on;
    plot(years(clusters == 1), precip_obs_1yr_stdnorm(clusters == 1),'or',...
        'LineStyle','none','MarkerFaceColor','r');
    plot(years(clusters == 2), precip_obs_1yr_stdnorm(clusters == 2),'ok',...
        'LineStyle','none','MarkerFaceColor','k');
    plot(years(clusters == 3), precip_obs_1yr_stdnorm(clusters == 3),'ob',...
        'LineStyle','none','MarkerFaceColor','b');
    ylm = get(gca,'Ylim');
    ylabel('Total Precipitation');
    
    subplot(7,5,20+10);
    y = linspace(ylm(1),ylm(2),200);
    plot(weights(1)*normpdf(y,mu_hat(1),sigma_hat(1)),y,'-r');
    hold on;
    plot(weights(2)*normpdf(y,mu_hat(2),sigma_hat(2)),y,'-k');
    plot(weights(3)*normpdf(y,mu_hat(3),sigma_hat(3)),y,'-b');
    ylim(ylm);

    % % % % % % % % % 
    % Now the occ data:
    [clusters, mu_hat, sigma_hat, weights] = GMM_cluster(occ_obs_1yr_stdnorm, occ_sim_1yr_stdnorm);
    subplot(7,5,20+(11:14));
    plot(years, occ_obs_1yr_stdnorm,'-k');
    hold on;
    plot(years(clusters == 1), occ_obs_1yr_stdnorm(clusters == 1),'or',...
        'LineStyle','none','MarkerFaceColor','r');
    plot(years(clusters == 2), occ_obs_1yr_stdnorm(clusters == 2),'ok',...
        'LineStyle','none','MarkerFaceColor','k');
    plot(years(clusters == 3), occ_obs_1yr_stdnorm(clusters == 3),'ob',...
        'LineStyle','none','MarkerFaceColor','b');
    ylm = get(gca,'Ylim');
    ylabel({'Precipitation','Occurrence'});
    
    subplot(7,5,20+15);
    y = linspace(ylm(1),ylm(2),200);
    plot(weights(1)*normpdf(y,mu_hat(1),sigma_hat(1)),y,'-r');
    hold on;
    plot(weights(2)*normpdf(y,mu_hat(2),sigma_hat(2)),y,'-k');
    plot(weights(3)*normpdf(y,mu_hat(3),sigma_hat(3)),y,'-b');
    ylim(ylm);
    %%
    
    
      
    outfilename = sprintf('plots/All_TS_%s_1yr_DOY%i_wet_season.png',stn_id,DOY);
    print(gcf,'-dpng',outfilename);
    
    %%
    % % %% % %% % % %% 
    % 2-yearly data:
    filename = sprintf('LL_%s_accum_DOY%i.mat',stn_id,DOY);
    load(filename);

    filename = sprintf('precip_%s_accum_DOY%i.mat',stn_id,DOY);
    load(filename);

    
    % Cluster the obs data:
    [clusters, mu_hat, sigma_hat, weights] = GMM_cluster(LL_obs_2yr_stdnorm, LL_sim_2yr_stdnorm);

    N = length(LL_obs_2yr_stdnorm);
    years = (2014.5-2*N):2:2012.5;
    if (DOY > 1) % using something later than Jan1 start date
        years = years+1; % so that it ends with 2014
    end
    
    figure;
    subplot(7,5,1:4);
    plot(years, -LL_obs_2yr_stdnorm,'-k');
    hold on;
    plot(years(clusters == 1), -LL_obs_2yr_stdnorm(clusters == 1),'ob',...
        'LineStyle','none','MarkerFaceColor','b');
    plot(years(clusters == 2), -LL_obs_2yr_stdnorm(clusters == 2),'ok',...
        'LineStyle','none','MarkerFaceColor','k');
    plot(years(clusters == 3), -LL_obs_2yr_stdnorm(clusters == 3),'or',...
        'LineStyle','none','MarkerFaceColor','r');
    ylm = get(gca,'Ylim');
    title(sprintf('Stn #%s: 2-yearly, Start DOY = %i',stn_id,DOY));
    ylabel('-Log-Likelihood');
    
    subplot(7,5,5);
    y = linspace(ylm(1),ylm(2),200);
    plot(weights(1)*normpdf(y,mu_hat(1),sigma_hat(1)),-y,'-b');
    hold on;
    plot(weights(2)*normpdf(y,mu_hat(2),sigma_hat(2)),-y,'-k');
    plot(weights(3)*normpdf(y,mu_hat(3),sigma_hat(3)),-y,'-r');
    ylim(ylm);
    
    % % % % % % % % % 
    % Now the precip data:
    [clusters, mu_hat, sigma_hat, weights] = GMM_cluster(precip_obs_2yr_stdnorm, precip_sim_2yr_stdnorm);
    subplot(7,5,6:9);
    plot(years, precip_obs_2yr_stdnorm,'-k');
    hold on;
    plot(years(clusters == 1), precip_obs_2yr_stdnorm(clusters == 1),'or',...
        'LineStyle','none','MarkerFaceColor','r');
    plot(years(clusters == 2), precip_obs_2yr_stdnorm(clusters == 2),'ok',...
        'LineStyle','none','MarkerFaceColor','k');
    plot(years(clusters == 3), precip_obs_2yr_stdnorm(clusters == 3),'ob',...
        'LineStyle','none','MarkerFaceColor','b');
    ylm = get(gca,'Ylim');
    ylabel('Total Precipitation');
    
    subplot(7,5,10);
    y = linspace(ylm(1),ylm(2),200);
    plot(weights(1)*normpdf(y,mu_hat(1),sigma_hat(1)),y,'-r');
    hold on;
    plot(weights(2)*normpdf(y,mu_hat(2),sigma_hat(2)),y,'-k');
    plot(weights(3)*normpdf(y,mu_hat(3),sigma_hat(3)),y,'-b');
    ylim(ylm);

    % % % % % % % % % 
    % Now the occ data:
    [clusters, mu_hat, sigma_hat, weights] = GMM_cluster(occ_obs_2yr_stdnorm, occ_sim_2yr_stdnorm);
    subplot(7,5,11:14);
    plot(years, occ_obs_2yr_stdnorm,'-k');
    hold on;
    plot(years(clusters == 1), occ_obs_2yr_stdnorm(clusters == 1),'or',...
        'LineStyle','none','MarkerFaceColor','r');
    plot(years(clusters == 2), occ_obs_2yr_stdnorm(clusters == 2),'ok',...
        'LineStyle','none','MarkerFaceColor','k');
    plot(years(clusters == 3), occ_obs_2yr_stdnorm(clusters == 3),'ob',...
        'LineStyle','none','MarkerFaceColor','b');
    ylm = get(gca,'Ylim');
    ylabel({'Precipitation','Occurrence'});
    
    subplot(7,5,15);
    y = linspace(ylm(1),ylm(2),200);
    plot(weights(1)*normpdf(y,mu_hat(1),sigma_hat(1)),y,'-r');
    hold on;
    plot(weights(2)*normpdf(y,mu_hat(2),sigma_hat(2)),y,'-k');
    plot(weights(3)*normpdf(y,mu_hat(3),sigma_hat(3)),y,'-b');
    ylim(ylm);
    
    %%
    filename = sprintf('LL_%s_accum_DOY%i_wet_season_only.mat',stn_id,DOY);
    load(filename);

    filename = sprintf('precip_%s_accum_DOY%i_wet_season_only.mat',stn_id,DOY);
    load(filename);

    
    % Cluster the obs data:
    [clusters, mu_hat, sigma_hat, weights] = GMM_cluster(LL_obs_2yr_stdnorm, LL_sim_2yr_stdnorm);

    N = length(LL_obs_2yr_stdnorm);
    years = (2014.5-2*N):2:2012.5;
    if (DOY > 1) % using something later than Jan1 start date
        years = years+1; % so that it ends with 2014
    end
    

    subplot(7,5,20+(1:4));
    plot(years, -LL_obs_2yr_stdnorm,'-k');
    hold on;
    plot(years(clusters == 1), -LL_obs_2yr_stdnorm(clusters == 1),'ob',...
        'LineStyle','none','MarkerFaceColor','b');
    plot(years(clusters == 2), -LL_obs_2yr_stdnorm(clusters == 2),'ok',...
        'LineStyle','none','MarkerFaceColor','k');
    plot(years(clusters == 3), -LL_obs_2yr_stdnorm(clusters == 3),'or',...
        'LineStyle','none','MarkerFaceColor','r');
    ylm = get(gca,'Ylim');
    title(sprintf('Stn #%s: 2-yearly, Start DOY = %i -- WET SEASON ONLY',stn_id,DOY));
    ylabel('-Log-Likelihood');
    
    subplot(7,5,20+5);
    y = linspace(ylm(1),ylm(2),200);
    plot(weights(1)*normpdf(y,mu_hat(1),sigma_hat(1)),-y,'-b');
    hold on;
    plot(weights(2)*normpdf(y,mu_hat(2),sigma_hat(2)),-y,'-k');
    plot(weights(3)*normpdf(y,mu_hat(3),sigma_hat(3)),-y,'-r');
    ylim(ylm);
    
    % % % % % % % % % 
    % Now the precip data:
    [clusters, mu_hat, sigma_hat, weights] = GMM_cluster(precip_obs_2yr_stdnorm, precip_sim_2yr_stdnorm);
    subplot(7,5,20+(6:9));
    plot(years, precip_obs_2yr_stdnorm,'-k');
    hold on;
    plot(years(clusters == 1), precip_obs_2yr_stdnorm(clusters == 1),'or',...
        'LineStyle','none','MarkerFaceColor','r');
    plot(years(clusters == 2), precip_obs_2yr_stdnorm(clusters == 2),'ok',...
        'LineStyle','none','MarkerFaceColor','k');
    plot(years(clusters == 3), precip_obs_2yr_stdnorm(clusters == 3),'ob',...
        'LineStyle','none','MarkerFaceColor','b');
    ylm = get(gca,'Ylim');
    ylabel('Total Precipitation');
    
    subplot(7,5,20+10);
    y = linspace(ylm(1),ylm(2),200);
    plot(weights(1)*normpdf(y,mu_hat(1),sigma_hat(1)),y,'-r');
    hold on;
    plot(weights(2)*normpdf(y,mu_hat(2),sigma_hat(2)),y,'-k');
    plot(weights(3)*normpdf(y,mu_hat(3),sigma_hat(3)),y,'-b');
    ylim(ylm);

    % % % % % % % % % 
    % Now the occ data:
    [clusters, mu_hat, sigma_hat, weights] = GMM_cluster(occ_obs_2yr_stdnorm, occ_sim_2yr_stdnorm);
    subplot(7,5,20+(11:14));
    plot(years, occ_obs_2yr_stdnorm,'-k');
    hold on;
    plot(years(clusters == 1), occ_obs_2yr_stdnorm(clusters == 1),'or',...
        'LineStyle','none','MarkerFaceColor','r');
    plot(years(clusters == 2), occ_obs_2yr_stdnorm(clusters == 2),'ok',...
        'LineStyle','none','MarkerFaceColor','k');
    plot(years(clusters == 3), occ_obs_2yr_stdnorm(clusters == 3),'ob',...
        'LineStyle','none','MarkerFaceColor','b');
    ylm = get(gca,'Ylim');
    ylabel({'Precipitation','Occurrence'});
    
    subplot(7,5,20+15);
    y = linspace(ylm(1),ylm(2),200);
    plot(weights(1)*normpdf(y,mu_hat(1),sigma_hat(1)),y,'-r');
    hold on;
    plot(weights(2)*normpdf(y,mu_hat(2),sigma_hat(2)),y,'-k');
    plot(weights(3)*normpdf(y,mu_hat(3),sigma_hat(3)),y,'-b');
    ylim(ylm);
    %%
    outfilename = sprintf('plots/All_TS_%s_2yr_DOY%i_wet_season.png',stn_id,DOY);
    print(gcf,'-dpng',outfilename);
    %%
    
    % % %% % %% % % %%
    % 3-yearly data:
    filename = sprintf('LL_%s_accum_DOY%i.mat',stn_id,DOY);
    load(filename);

    filename = sprintf('precip_%s_accum_DOY%i.mat',stn_id,DOY);
    load(filename);

    % Cluster the obs data:
    [clusters, mu_hat, sigma_hat, weights] = GMM_cluster(LL_obs_3yr_stdnorm, LL_sim_3yr_stdnorm);

    N = length(LL_obs_3yr_stdnorm);
    years = (2015-3*N):3:2012;
    if (DOY > 1) % using something later than Jan1 start date
        years = years+1; % so that it ends with 2014
    end
    
    figure;
    subplot(7,5,1:4);
    plot(years, -LL_obs_3yr_stdnorm,'-k');
    hold on;
    plot(years(clusters == 1), -LL_obs_3yr_stdnorm(clusters == 1),'ob',...
        'LineStyle','none','MarkerFaceColor','b');
    plot(years(clusters == 2), -LL_obs_3yr_stdnorm(clusters == 2),'ok',...
        'LineStyle','none','MarkerFaceColor','k');
    plot(years(clusters == 3), -LL_obs_3yr_stdnorm(clusters == 3),'or',...
        'LineStyle','none','MarkerFaceColor','r');
    ylm = get(gca,'Ylim');
    title(sprintf('Stn #%s: 3-yearly, Start DOY = %i',stn_id,DOY));
    ylabel('-Log-Likelihood');
   
    subplot(7,5,5);
    y = linspace(ylm(1),ylm(2),200);
    plot(weights(1)*normpdf(y,mu_hat(1),sigma_hat(1)),-y,'-b');
    hold on;
    plot(weights(2)*normpdf(y,mu_hat(2),sigma_hat(2)),-y,'-k');
    plot(weights(3)*normpdf(y,mu_hat(3),sigma_hat(3)),-y,'-r');
    ylim(ylm);

    % % % % % % % % % 
    % Now the precip data:
    [clusters, mu_hat, sigma_hat, weights] = GMM_cluster(precip_obs_3yr_stdnorm, precip_sim_3yr_stdnorm);
    subplot(7,5,6:9);
    plot(years, precip_obs_3yr_stdnorm,'-k');
    hold on;
    plot(years(clusters == 1), precip_obs_3yr_stdnorm(clusters == 1),'or',...
        'LineStyle','none','MarkerFaceColor','r');
    plot(years(clusters == 2), precip_obs_3yr_stdnorm(clusters == 2),'ok',...
        'LineStyle','none','MarkerFaceColor','k');
    plot(years(clusters == 3), precip_obs_3yr_stdnorm(clusters == 3),'ob',...
        'LineStyle','none','MarkerFaceColor','b');
    ylm = get(gca,'Ylim');
    ylabel('Total Precipitation');
    
    subplot(7,5,10);
    y = linspace(ylm(1),ylm(2),200);
    plot(weights(1)*normpdf(y,mu_hat(1),sigma_hat(1)),y,'-r');
    hold on;
    plot(weights(2)*normpdf(y,mu_hat(2),sigma_hat(2)),y,'-k');
    plot(weights(3)*normpdf(y,mu_hat(3),sigma_hat(3)),y,'-b');
    ylim(ylm);

    % % % % % % % % % 
    % Now the occ data:
    [clusters, mu_hat, sigma_hat, weights] = GMM_cluster(occ_obs_3yr_stdnorm, occ_sim_3yr_stdnorm);
    subplot(7,5,11:14);
    plot(years, occ_obs_3yr_stdnorm,'-k');
    hold on;
    plot(years(clusters == 1), occ_obs_3yr_stdnorm(clusters == 1),'or',...
        'LineStyle','none','MarkerFaceColor','r');
    plot(years(clusters == 2), occ_obs_3yr_stdnorm(clusters == 2),'ok',...
        'LineStyle','none','MarkerFaceColor','k');
    plot(years(clusters == 3), occ_obs_3yr_stdnorm(clusters == 3),'ob',...
        'LineStyle','none','MarkerFaceColor','b');
    ylm = get(gca,'Ylim');
    ylabel({'Precipitation','Occurrence'});
    
    subplot(7,5,15);
    y = linspace(ylm(1),ylm(2),200);
    plot(weights(1)*normpdf(y,mu_hat(1),sigma_hat(1)),y,'-r');
    hold on;
    plot(weights(2)*normpdf(y,mu_hat(2),sigma_hat(2)),y,'-k');
    plot(weights(3)*normpdf(y,mu_hat(3),sigma_hat(3)),y,'-b');
    ylim(ylm);
    
    %%
    filename = sprintf('LL_%s_accum_DOY%i_wet_season_only.mat',stn_id,DOY);
    load(filename);

    filename = sprintf('precip_%s_accum_DOY%i_wet_season_only.mat',stn_id,DOY);
    load(filename);
    % Cluster the obs data:
    [clusters, mu_hat, sigma_hat, weights] = GMM_cluster(LL_obs_3yr_stdnorm, LL_sim_3yr_stdnorm);

    N = length(LL_obs_3yr_stdnorm);
    years = (2015-3*N):3:2012;
    if (DOY > 1) % using something later than Jan1 start date
        years = years+1; % so that it ends with 2014
    end
    
    subplot(7,5,20+(1:4));
    plot(years, -LL_obs_3yr_stdnorm,'-k');
    hold on;
    plot(years(clusters == 1), -LL_obs_3yr_stdnorm(clusters == 1),'ob',...
        'LineStyle','none','MarkerFaceColor','b');
    plot(years(clusters == 2), -LL_obs_3yr_stdnorm(clusters == 2),'ok',...
        'LineStyle','none','MarkerFaceColor','k');
    plot(years(clusters == 3), -LL_obs_3yr_stdnorm(clusters == 3),'or',...
        'LineStyle','none','MarkerFaceColor','r');
    ylm = get(gca,'Ylim');
    title(sprintf('Stn #%s: 3-yearly, Start DOY = %i -- WET SEASON ONLY',stn_id,DOY));
    ylabel('-Log-Likelihood');
   
    subplot(7,5,20+5);
    y = linspace(ylm(1),ylm(2),200);
    plot(weights(1)*normpdf(y,mu_hat(1),sigma_hat(1)),-y,'-b');
    hold on;
    plot(weights(2)*normpdf(y,mu_hat(2),sigma_hat(2)),-y,'-k');
    plot(weights(3)*normpdf(y,mu_hat(3),sigma_hat(3)),-y,'-r');
    ylim(ylm);

    % % % % % % % % % 
    % Now the precip data:
    [clusters, mu_hat, sigma_hat, weights] = GMM_cluster(precip_obs_3yr_stdnorm, precip_sim_3yr_stdnorm);
    subplot(7,5,20+(6:9));
    plot(years, precip_obs_3yr_stdnorm,'-k');
    hold on;
    plot(years(clusters == 1), precip_obs_3yr_stdnorm(clusters == 1),'or',...
        'LineStyle','none','MarkerFaceColor','r');
    plot(years(clusters == 2), precip_obs_3yr_stdnorm(clusters == 2),'ok',...
        'LineStyle','none','MarkerFaceColor','k');
    plot(years(clusters == 3), precip_obs_3yr_stdnorm(clusters == 3),'ob',...
        'LineStyle','none','MarkerFaceColor','b');
    ylm = get(gca,'Ylim');
    ylabel('Total Precipitation');
    
    subplot(7,5,20+10);
    y = linspace(ylm(1),ylm(2),200);
    plot(weights(1)*normpdf(y,mu_hat(1),sigma_hat(1)),y,'-r');
    hold on;
    plot(weights(2)*normpdf(y,mu_hat(2),sigma_hat(2)),y,'-k');
    plot(weights(3)*normpdf(y,mu_hat(3),sigma_hat(3)),y,'-b');
    ylim(ylm);

    % % % % % % % % % 
    % Now the occ data:
    [clusters, mu_hat, sigma_hat, weights] = GMM_cluster(occ_obs_3yr_stdnorm, occ_sim_3yr_stdnorm);
    subplot(7,5,20+(11:14));
    plot(years, occ_obs_3yr_stdnorm,'-k');
    hold on;
    plot(years(clusters == 1), occ_obs_3yr_stdnorm(clusters == 1),'or',...
        'LineStyle','none','MarkerFaceColor','r');
    plot(years(clusters == 2), occ_obs_3yr_stdnorm(clusters == 2),'ok',...
        'LineStyle','none','MarkerFaceColor','k');
    plot(years(clusters == 3), occ_obs_3yr_stdnorm(clusters == 3),'ob',...
        'LineStyle','none','MarkerFaceColor','b');
    ylm = get(gca,'Ylim');
    ylabel({'Precipitation','Occurrence'});
    
    subplot(7,5,20+15);
    y = linspace(ylm(1),ylm(2),200);
    plot(weights(1)*normpdf(y,mu_hat(1),sigma_hat(1)),y,'-r');
    hold on;
    plot(weights(2)*normpdf(y,mu_hat(2),sigma_hat(2)),y,'-k');
    plot(weights(3)*normpdf(y,mu_hat(3),sigma_hat(3)),y,'-b');
    ylim(ylm);
    %%
    
    outfilename = sprintf('plots/All_TS_%s_3yr_DOY%i_wet_season.png',stn_id,DOY);
    print(gcf,'-dpng',outfilename);
    
end % for loop over stations
    
end % function
