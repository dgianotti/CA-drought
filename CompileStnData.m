
clear;
clc;

if strcmpi(getenv('OS'),'Windows_NT')
    addpath('C:\Users\gianotti\Documents\IntensityLib\')
else % linux
    addpath('IntensityLib');
end

today_doy = today - datenum(year(today),1,1) + 1;
remove_seasonal_cycle = false;

%% Download and calculate all of the LL data:
calculate_daily_LL_data;

%% Create different accumulated LL time series (each gets saved with the 
% daily LL data for that stn).

% Delete old accumulated data files:
delete('*accum_DOY*.mat');
accumulate_LL({'oneyear','twoyear','threeyear'},today_doy,remove_seasonal_cycle);
accumulate_LL({'oneyear','twoyear','threeyear'},1,remove_seasonal_cycle);

% Transform the annualized (or bi-annualized, etc.) data to be more normal 
% using the sim distribution
cdf_transform_accum_LL_data();

% Cluster the obs data and make some plots:
make_accumlated_cluster_plots_LL(1);
make_accumlated_cluster_plots_LL(today_doy);

%% Now do the same thing for total precipitation and occurrence instead of LL:

accumulate_precip({'oneyear','twoyear','threeyear'},today_doy,remove_seasonal_cycle);
accumulate_precip({'oneyear','twoyear','threeyear'},1,remove_seasonal_cycle);

% Transform the annualized (or bi-annualized, etc.) data to be more normal 
% using the sim distribution
cdf_transform_accum_precip_data();

% Cluster the obs data and make some plots:
make_accumlated_cluster_plots(1);
make_accumlated_cluster_plots(today_doy);

% Make annual LL/precip scatter plots and seasonal time-series:
make_annual_LL_precip_plots;

%% Try a comparison between May 1 start date with and without May-Oct precip:
may1 = 121;
oct31 = 304;
accumulate_LL({'oneyear','twoyear','threeyear'},may1,0);
accumulate_precip({'oneyear','twoyear','threeyear'},may1,0);

accumulate_precip_wet_season_only({'oneyear','twoyear','threeyear'},may1,0);
accumulate_LL_wet_season_only({'oneyear','twoyear','threeyear'},may1,0);

cdf_transform_accum_LL_data;
cdf_transform_accum_precip_data;

make_accumulated_cluster_plots_wet_season(may1);


%%

% Convert the sim to a normal distribution using order statistics:
[std_norm_pdfs_sim, uniform_pdfs_sim, ~] = empirical_2_normal_via_order_stats(LL_sim_annual);

% Get the std_normal values for the obs using the sim transform as a
% look-up table:
std_norm_pdfs_obs = sim_dist_to_obs_std_norm(LL_obs_annual, LL_sim_annual, ...
    uniform_pdfs_sim);


subplot(3,5,i);
plot(years,quantile(LL_sim_annual,[.025,.5,.975],2),'-r');
hold on;
plot(years,LL_obs_annual,'-b');
xlim([10*floor(min(years/10)),2020]);
title(id);
set(gca,'XTick',1900:50:2000);


% subplot(3,5,i);
% plot(years,quantile(std_norm_pdfs_sim,[.025,.5,.975],2),'-r');
% hold on;
% plot(years,std_norm_pdfs_obs,'-b');
% xlim([10*floor(min(years/10)),2020]);
% title(id);
% set(gca,'XTick',[1900:50:2000]);
end
print(gcf,'-dpdf','DailyNorm_AnnualNotNorm_Apr2Apr.pdf');

%% Make annual plots!
    
for i = 1:length(good_CA_IDs)
    id = good_CA_IDs{i};
    load(['LL_',id,'.mat']);
    
    bad_data_obs = isnan(LL_obs);
    LL_obs_annual = nanmean(LL_obs,2);

    n_years = length(years);
    num_useable_sims = floor(length(LL_sim)/n_years);
    bad_data_sim = repmat(bad_data_obs,[num_useable_sims,1]);
       
    % Put NaNs into the sim data at the same locations...
    LL_sim = LL_sim(1:(n_years*num_useable_sims),:);
    LL_sim(bad_data_sim) = nan;
    LL_sim = nanmean(LL_sim,2);
    
    LL_sim = reshape(LL_sim, [n_years,num_useable_sims]);
    
    plot(years,LL_obs_annual,'-k');
    hold on;
    plot(years,quantile(LL_sim,.05,2),'-b');
    plot(years,quantile(LL_sim,.95,2),'-b');
    plot(years,quantile(LL_sim,.5,2),'-r');
    ylabel('Log Likelihood');
    title(sprintf('Normalized annual log likelihood precipitation for station %s',id)); 
    
    print(gcf,'-dpng',sprintf('LL_annual_%s.png',id));
    close all;
end
%% Plot LL versus precip for each station:
% figure;
% 
% for i = 1:length(good_CA_IDs)
%     id = good_CA_IDs{i}
%     load(['LL_',id,'.mat']);
%     SimStn = load_stn_data(id,'SimStn');
%     ImpStn = load_stn_data(id,'ImpStn');
%     subplot(5,3,i);
%     sim_precip = SimStn.intensity_data(1:size(LL_sim,1),:);
%     obs_precip = ImpStn.intensity_data; % Ignoring new (post 2010) data for now...
%     LL_obs = LL_obs(1:size(obs_precip,1),:);
%     
%     scatter(sim_precip(:),LL_sim(:),'.r');
%     hold on;
%     scatter(obs_precip(:),LL_obs(:),'.b');        
%     title(id)
% end
% 
% print(gcf,'-dpng','LL_vs_precip.png');




%% Let's do some GMM fitting!

for i = 1:length(good_CA_IDs)
    close all;
    id = good_CA_IDs{i}
    filename = ['LL_',id,'.mat'];    
    load(filename); % loads LL_obs [N x 1], LL_sim [N x ~1000], years [N x 1]
    
    
    m1=mean(mean(LL_sim));
    sig1=mean(std(LL_sim'));
    poss=[m1-5*sig1:1:m1+5*sig1];

    % first try for sim data where should get one cluster only

    options = statset('MaxIter',1000,'Display','Final');
    replicates=10;
    
    for simulation=1:50
        simulation
        X = LL_sim(:,simulation);

        for num_clusters=1:3
            %gm = gmdistribution.fit(X,nnn,'regularize',[0.01],'Options',options,'CovType','full','Replicates',replicates);
            gm = gmdistribution.fit(X,num_clusters,'regularize',[0.01],'Options',options,'CovType','diag','Replicates',replicates);

            aicsave_sim(num_clusters,simulation)=gm.AIC;
            bicsave_sim(num_clusters,simulation)=gm.BIC;

        end

    end
    
    [~,bchoiceaic]=min(aicsave_sim);
    [~,bchoicebic]=min(bicsave_sim);
    figure

    hist([bchoiceaic',bchoicebic'])
    title('number of times each order chosen when applied to null (simulated likelihoods)')
    legend('AIC','BIC')

    X = LL_obs;

    replicates=50;

    for num_clusters=1:3
        %gm = gmdistribution.fit(X,nnn,'regularize',[0.01],'Options',options,'CovType','full','Replicates',replicates);
        gm = gmdistribution.fit(X,num_clusters,'regularize',[0.01],'Options',options,'CovType','diag','Replicates',replicates);

        aicsave(num_clusters,1)=gm.AIC;
        bicsave(num_clusters,1)=gm.BIC;

    end

    [~,min_BIC_idx_obs]=min(bicsave);
    selected_num_clusters = min_BIC_idx_obs

    figure
    plot(bicsave)
    hold on
    plot(min_BIC_idx_obs,bicsave(min_BIC_idx_obs),'ko','linewidth',[3])
    title('BIC for observed data showing optimal number of clusters')
    ylabel('BIC')
    xlabel('number of clusters')


    % redo fit for chosen number of clusters

    %gm = gmdistribution.fit(X,min_BIC_idx_obs,'Options',options,'CovType','full','Replicates',replicates);
    gm = gmdistribution.fit(X,min_BIC_idx_obs,'regularize',[0.01],'Options',options,'CovType','diag','Replicates',replicates);

    gm.NlogL
    
    gm.PComponents
    gm.mu
    gm.Sigma

    pdfmat = zeros(length(poss),min_BIC_idx_obs);
    for comp=1:min_BIC_idx_obs
        pdfmat(:,comp)=gm.PComponents(comp).*normpdf(poss,gm.mu(comp),gm.Sigma(comp).^0.5);
    end

    figure

    plot(poss,pdfmat)
    hold on
    title('component pdfs')
    xlabel('likelihood ')
    
    idx=cluster(gm,X);


    figure
    plot(X)
    hold on
    II=find(idx==1);
    plot(II,X(II),'ro','linewidth',[3]);
    
    II=find(idx==2);
    plot(II,X(II),'go','linewidth',[3]);
    
    II=find(idx==3);
    plot(II,X(II),'bo','linewidth',[3]);
    
    II=find(idx==4);
    plot(II,X(II),'ko','linewidth',[3]);
    
    II=find(idx==5);
    plot(II,X(II),'mo','linewidth',[3]);
    
    legend('data','cluster 1','cluster 2')


    % home made version
    m1=mean(mean(LL_sim));
    sig1=mean(std(LL_sim'));
    poss=[m1-5*sig1:1:m1+5*sig1];
    
    X = LL_obs;
    
    %mu1 = median(X);
    mu1 = mean(X);
    
    sig1=mean(std(LL_sim'));
    
    
    lowendmu=(mu1-4*sig1):(3*sig1/10):(mu1-1*sig1);
    highendmu=(mu1+1*sig1):(3*sig1/10):(mu1+4*sig1);
    
    lowendsig=sig1.*(0.1:0.025:0.2);
    highendsig=sig1.*(0.1:0.025:0.2);
    
    mix1s = 0.8:0.025:1;
    mix2s = 0:0.025:0.1;
        
    cc=0;
    numlen=length(lowendmu).*length(lowendsig).*length(highendmu).*length(highendsig).*length(mix1s).*length(mix2s)
    
    parsave=zeros(numlen,9);
    llsave=zeros(numlen,1);
    pensave=zeros(numlen,1);

    for mu2=lowendmu
        for sig2=lowendsig
            for mu3=highendmu
                for sig3=highendsig
                    for mix1=mix1s
                        
                        for mix2=mix2s
                            if mix2+mix1<=1
                                cc=cc+1;
                                
                                percentdone=round(100*(cc./numlen));
                                
                                if mod(cc,1000)==0
                                    percentdone
                                end
                                
                                
                                mix3=1-(mix1+mix2);
                                ll=log(mix1.*normpdf(X,mu1,sig1)+mix2.*normpdf(X,mu2,sig2)+mix3.*normpdf(X,mu3,sig3));
                                
                                llsave(cc,1)=sum(ll);
                                parsave(cc,:)=[mu1,sig1,mu2,sig2,mu3,sig3,mix1,mix2,mix3];
                                
                                pensave(cc,1)=(mix1~=0).*3+(mix2~=0).*3+(mix3~=0).*3;
                                
                                
                            end
                            
                        end
                    end
                end
            end
        end
    end

    [~,min_BIC_idx_obs]=sort(-2*llsave(1:cc)+2*pensave(1:cc));
    

    II=find((parsave(1:cc,8)==0)&(parsave(1:cc,9)>0));
    
    [aa1,bb1]=min(-2*llsave(II)+2*pensave(II));
    bb1A=II(bb1)
    
    
    II=find((parsave(1:cc,9)==0)&(parsave(1:cc,8)>0));
    
    [aa2,bb2]=min(-2*llsave(II)+2*pensave(II));
    bb2A=II(bb2)
    
    
    II=find((parsave(1:cc,8)==0)&(parsave(1:cc,9)==0));
    
    [aa3,bb3]=min(-2*llsave(II)+2*pensave(II));
    
    bb3A=II(bb3)
    
    [aaa,bbb]=min([aa1,aa2,aa3])
    
    bvals=[bb1A,bb2A,bb3A]
    
    b=bvals(bbb)
    
    
    
    
    parsave(b,:)
    mu1=parsave(b,1);
    sig1=parsave(b,2);
    mu2=parsave(b,3);
    sig2=parsave(b,4);
    mu3=parsave(b,5);
    sig3=parsave(b,6);
    mix1=parsave(b,7);
    mix2=parsave(b,8);
    mix3=parsave(b,9);
    
    [mu1,sig1,mu2,sig2,mu3,sig3,mix1,mix2,mix3]
    
    poss=[mu1-5*sig1:1:mu1+5*sig1];
    
    pdf1=mix1.*normpdf(poss,mu1,sig1);
    pdf2=mix2.*normpdf(poss,mu2,sig2);
    pdf3=+mix3.*normpdf(poss,mu3,sig3);
    pdfall=mix1.*normpdf(poss,mu1,sig1)+mix2.*normpdf(poss,mu2,sig2)+mix3.*normpdf(poss,mu3,sig3);
    
    pdfmat(:,1)=pdf1;
    pdfmat(:,2)=pdf2;
    pdfmat(:,3)=pdf3;
    
    
    
    figure
    
    plot(poss,pdfmat)
    hold on
    title('component pdfs')
    xlabel('likelihood ')
    
    
    
    plot(poss,pdf1,'g','linewidth',3)
    hold on
    plot(poss,pdf2,'r','linewidth',3)
    plot(poss,pdf3,'b','linewidth',3)
    plot(poss,pdfall,'k','linewidth',3)
    
    


    % now for the associateion part...
    
    
    pdf1dat=mix1.*normpdf(X,mu1,sig1);
    pdf2dat=mix2.*normpdf(X,mu2,sig2);
    pdf3dat=+mix3.*normpdf(X,mu3,sig3);
    
    [aa,idx]=max([pdf1dat,pdf2dat,pdf3dat]');
    
    % idx has in it the cluster of each data point
    
    figure
    plot(X)
    hold on
    II = find(idx==1);
    plot(II,X(II),'ro','linewidth',[3]);
    
    II=find(idx==2);
    plot(II,X(II),'go','linewidth',[3]);
    
    II=find(idx==3);
    plot(II,X(II),'bo','linewidth',[3]);
    
    II=find(idx==4);
    plot(II,X(II),'ko','linewidth',[3]);
    
    II=find(idx==5);
    plot(II,X(II),'mo','linewidth',[3]);
    
    if sum(idx==2)>0
        legend('data','cluster 1','cluster 2')
    end
    
    if sum(idx==2)>0 && sum(idx==3)>0
        legend('data','cluster 1','cluster 2','cluster 3')
    end
    
    
    keyboard;
end