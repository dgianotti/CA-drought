% Let's get the (current) data for the CA stations!

clear;
clc;

addpath('C:\Users\gianotti\Documents\IntensityLib');

load('CA_ids.mat');

fraction_missing_vec = [ 0.2302, 0.3140, 1, 0.2133, 0.4727, ...
    0.0111, 0.0364, 0.0078, 0.2055, 0.0559, ...
    0.6050, 0.0013, 0.0124, 0.3296, 0.2581, ...
    0.1248, 0.0377, 0.0013, 0.2302, 0.0663, ...
    0.0026, 0.0020, 0.0104, 0.0228, 0.0449, ...
    0.1990, 0.5379, 0.1586, 0.1944, 0.5879, ...
    1, 0.0845, 0.1118, 0.1190, 0.0059 ];

good_CA_IDs = CA_IDs(fraction_missing_vec < 0.05);
    

%% Loop over stations:
for i = 1:length(good_CA_IDs)

    id = good_CA_IDs{i}
    % Download the latest USHCN data:
    url = ['http://www1.ncdc.noaa.gov/pub/data/ghcn/daily/hcn/USC00',id,'.dly'];
    filename = ['GHCN-Daily/',id,'.dly'];
    urlwrite(url,filename);
    
    
    % Read in fixed-width data:
    [datenums,precip] = format_GHCN_precip_data(filename,'PadFirstLastYears','ExcludeLeapDays');

    %     % Set last day to to:
    %     feb28_2014 = datenum('2014-02-28','yyyy-mm-dd');
    %     datenums = [datenums(:); ((datenums(end)+1):feb28_2014)'];
    %     precip = [precip(:); nan([length(datenums)-length(precip),1])];
        
    new_precip = precip(datenums>= datenum('2010-01-01') & datenums<= today);
    new_datenums = datenums(datenums>= datenum('2010-01-01') & datenums<= today);

    % Determine fraction of data missing since 2010:
    fraction_missing = sum(isnan(new_precip))/numel(new_precip)

    
    %if (mod(length(new_precip),365) ~= 0) %uh oh... it should...
    %error('The new data should be a multiple of 365 days... did you use ExcludeLeapDays in format_GHCN_precip_data\n');
    %end


    % Fill missing data values somehow!
    % Neighboring GHCN sites?
    missing_data_datenums = new_datenums(isnan(new_precip)) ;
    
    % Look up lat/lon in ushcn-stations.txt:
    [lat, lon] = get_ushcn_lat_lon(id);
    
    % Maybe use GSOD:
    % ftp://ftp.ncdc.noaa.gov/pub/data/gsod
    
    % filled_precip = get_precip_GSOD(missing_data_datenums,lat,lon,20);

    % Or NOAA's Quality controlled local climatological data:
    % http://cdo.ncdc.noaa.gov/qclcd/QCLCD
    % filled_precip = get_precip_QCLCD(missing_data_datenums,lat,lon);

    % Load old precip data:
    ImpStn = load_stn_data(id,'ImpStn');

    % Determine likelihood of data given model:
    
    nan_padded_new_data = nan(365*5,1);
    nan_padded_new_data(1:length(new_precip)) = new_precip;
    
    data = [ImpStn.intensity_data; reshape(nan_padded_new_data,[365,5])'];
        
    normalize_LL = true;
    start_year = 2010 - ImpStn.num_years;
    %[LL_obs, years] = get_annual_log_likelihood(id, data, start_year, normalize_LL);
    [LL_obs, years] = get_daily_log_likelihood(id, data, start_year);
    
    % Now sim data!
    SimStn = load_stn_data(id,'SimStn');        

    %[LL_sim, ~] = get_annual_log_likelihood(id, SimStn.intensity_data, start_year, normalize_LL);
    [LL_sim, ~] = get_daily_log_likelihood(id, SimStn.intensity_data, start_year);
    save(sprintf('LL_%s.mat',id),'LL_obs','LL_sim','years');
    
end 
    

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


%% Plot annual cycle of precip and LL at each station:


for i = 1:length(good_CA_IDs)
    figure;
    id = good_CA_IDs{i}
    load(['LL_',id,'.mat']);
    SimStn = load_stn_data(id,'SimStn');
    ImpStn = load_stn_data(id,'ImpStn');
    subplot(2,2,1);
    sim_precip = SimStn.intensity_data(1:size(LL_sim,1),:);
   
    obs_precip = ImpStn.intensity_data; % Ignoring new (post 2010) data for now...
    LL_obs = LL_obs(1:size(obs_precip,1),:);
    
    plot(1:365,mean(sim_precip,1),'-r');
    hold on;
    plot(1:365,mean(obs_precip,1),'-b');
    xlabel('Day of Year');
    xlim([1,365]);
    ylabel('Mean Daily Precip. [mm]');
    legend({'Sim','Obs'});
    title(id);
    
    subplot(2,2,3);
    plot(1:365,mean(LL_sim,1),'-r');
    hold on;
    plot(1:365,nanmean(LL_obs,1),'-b');
    xlabel('Day of Year');
    xlim([1,365]);
    ylabel('Mean Daily LL');
    legend({'Sim','Obs'});

    subplot(2,2,[2,4]);
    sim_precip = SimStn.intensity_data(1:size(LL_sim,1),:);
    obs_precip = ImpStn.intensity_data(1:size(LL_obs,1),:);
    scatter(sim_precip(:),LL_sim(:),'.r');
    hold on;
    scatter(obs_precip(:),LL_obs(:),'.b');        
    xlabel('Precipitation [mm]');
    ylabel('Log-Likelihood');
    legend({'Sim','Obs'});
    
    print(gcf,'-dpng',sprintf('Seasonal_precip_LL_%s.png',id));
    close all;
    
end





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