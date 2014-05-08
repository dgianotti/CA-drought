function make_annual_LL_precip_plots()

% Plot annual cycle of precip and LL at each station:

load('CA_ids.mat');

for i = 1:length(good_CA_IDs)
    figure;
    fprintf('Plotting figure for station %i of %i...\n',i,length(good_CA_IDs));
    id = good_CA_IDs{i};
    load(['LL_',id,'.mat']);
    
    SimStn = load_stn_data(id,'SimStn');
    
    subplot(6,2,[1,3]);
    precip_sim = SimStn.intensity_data;

    precip_obs(1,:) = [];
        
    plot(1:365,nanmean(precip_sim,1),'-r','LineWidth',2);
    hold on;
    plot(1:365,nanmean(precip_obs,1),'-b','LineWidth',2);
    xlabel('Day of Year');
    xlim([1,365]);
    ylabel('Mean Daily Precip. [mm]');
    legend({'Sim','Obs'});
    title(id);
    
    subplot(6,2,[5,7]);
    plot(1:365,nanmean(LL_sim,1),'-r','LineWidth',2);
    hold on;
    plot(1:365,nanmean(LL_obs,1),'-b','LineWidth',2);
    xlabel('Day of Year');
    xlim([1,365]);
    ylabel('Mean Daily LL');
    legend({'Sim','Obs'});

    subplot(6,2,[9,11]);
    plot(1:365,nanmean(precip_sim>0,1),'-r','LineWidth',2);
    hold on;
    plot(1:365,nanmean(precip_obs>0,1),'-b','LineWidth',2);
    xlabel('Day of Year');
    xlim([1,365]);
    ylabel({'Mean Daily', 'Occurrence Frequency'});
    legend({'Sim','Obs'});

    
    subplot(6,2,[2,4,6]);
    precip_sim = precip_sim(1:size(LL_sim,1),:);
    plot(precip_sim(:),LL_sim(:),'.r','LineStyle','none');
    hold on;
    plot(precip_obs(:),LL_obs(:),'.b','LineStyle','none');        
    xlabel('Precipitation [mm]');
    ylabel('Log-Likelihood');
    legend({'Sim','Obs'});
    
    subplot(6,2,[8,10,12]);
    plot(precip_sim(:)>0,LL_sim(:),'.r','LineStyle','none');
    hold on;
    plot(precip_obs(:)>0,LL_obs(:),'.b','LineStyle','none');        
    xlabel('Precipitation [mm]');
    ylabel('Log-Likelihood');
    legend({'Sim','Obs'});

    print(gcf,'-dpng',sprintf('plots/Seasonal_precip_LL_%s.png',id));
    close all;
    
end

