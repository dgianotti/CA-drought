% From Guido's home made version

clear;
clc;

load('CA_ids');

% Load the likelihood data:
stn = 1;
id = good_CA_IDs{stn};

load(['LL_',id,'.mat']);


% Get the annual likelihoods:
LL_sim = 365*nanmean(LL_sim,2)';
LL_obs = 365*nanmean(LL_obs,2)';

num_years = length(LL_obs);
num_sims = floor( length(LL_sim)/num_years );

% Round sim to a multiple of obs length
LL_sim( (length(LL_obs)*num_sims + 1):end ) = [];

%%


prob_obs = ksdensity(LL_sim,LL_obs,'function','cdf');

LL_obs_stdnorm = norminv(prob_obs,0,1);

mu = 0;
stdev = 1;
FORCEMEANSTD = 1;



LCRIT=0.001;

for num_clusters = 1:5
    disp(num_clusters);

    for ens = 1:100
        % randomly partition data into clusters
        assigned_cluster = zeros(num_years,num_sims);
        assigned_cluster(:,1) = ceil(rand(num_years,1)*num_clusters);

        % estimate parameters assuming ident correct
        plotit = 0;
        done = 0;

        for sim = 1:num_sims
            
            if done==0
                
                for clust = 1:num_clusters
                    % If its the first time through:
                    if sim == 1
                        % Find the years currently assigned to this cluster
                        II = find(assigned_cluster(:,sim)==clust);
                        mixes(clust,sim) = length(II)./num_years;
                    else
                        pr=pdfdat(:,clust)./sum(pdfdat,2);
                        mixes(clust,sim)=mean(pr);
                        means(clust,sim)=sum(pr.*LL_obs_stdnorm')./sum(pr);
                        %
                        %                     if ii==2
                        %                         means(ii,rep)=min(means(ii,rep),-1); % at least
                        %                         one std to left
                        %                     end
                        %
                        %                     if ii==3
                        %                         means(ii,rep)=max(means(ii,rep),1); % at least
                        %                         one std to left
                        %                     end
                        %
                        stds(clust,sim)=sqrt(sum(pr.*(LL_obs_stdnorm'-means(clust,sim)).^2)./sum(pr));
                        stds(clust,sim)=max(stds(clust,sim),0.1);
                        
                        if FORCEMEANSTD==1
                            if clust==1
                                means(clust,sim)=mu;
                                stds(clust,sim)=stdev;
                            end
                        end
                        
                    end
                    
                    % hard clust if rep==1
                    if sim==1
                        means(clust,sim)=mean(LL_obs_stdnorm(II));
                        stds(clust,sim)=std(LL_obs_stdnorm(II));
                    end
                end % loop over clusters
                
                % create the prob model, and re-assign ident
                
                pdfdat=zeros(num_years,num_clusters);
                for clust=1:num_clusters
                    
                    pdfdat(:,clust)=mixes(clust,sim).*normpdf(LL_obs_stdnorm,means(clust,sim),stds(clust,sim));
                end
                
                [aa,idxx]=max(pdfdat');
                
                assigned_cluster(:,sim+1)=idxx;
                
                if plotit==1
                    poss=[min(LL_obs_stdnorm)-2:0.1:max(LL_obs_stdnorm)+2];
                    
                    for clust=1:num_clusters
                        pdfmat(:,clust)=mixes(clust,sim).*normpdf(poss,means(clust,sim),stds(clust,sim));
                    end
                    
                    plot(poss,pdfmat(:,1),'g')
                    hold on
                    plot(poss,pdfmat(:,2),'r')
                    
                    plot(poss,sum(pdfmat,2),'k')
                    ksdensity(LL_obs_stdnorm)
                    
                end
                
                llsave(sim,1)=sum(log(sum(pdfdat,2)));
                
                if sim>1
                    if abs(llsave(sim,1)-llsave(sim-1,1))<LCRIT
                        done=1;
                        replast=sim;
                    end
                end
                
                
            end
            
        end
               
        AICs(num_clusters,ens)=-2*llsave(replast,1)+2*num_clusters*2+2*(num_clusters-1);
        loglikesave(num_clusters,ens)=-2*llsave(replast,1);
        repsave(num_clusters,ens)=replast;
        
        parsave(num_clusters,1:num_clusters,1,ens)=means(:,replast)';
        parsave(num_clusters,1:num_clusters,2,ens)=stds(:,replast)';
        parsave(num_clusters,1:num_clusters,3,ens)=mixes(:,replast)';
        
    end
end

for num_clusters=1:5
    [a,b]=min(AICs(num_clusters,:));
    AICs(num_clusters,1)=a;
    loglikesave(num_clusters,1)=loglikesave(num_clusters,b);
    parsave(num_clusters,1:num_clusters,1,1)=parsave(num_clusters,1:num_clusters,1,b);
    parsave(num_clusters,1:num_clusters,2,1)=parsave(num_clusters,1:num_clusters,2,b);
    parsave(num_clusters,1:num_clusters,3,1)=parsave(num_clusters,1:num_clusters,3,b);
    
    
end



[~,selected_num_clusters]=min(AICs(:,1));

MEANS = parsave(b,1:b,1,1)
STDS = parsave(b,1:b,2,1)
MIXES = parsave(b,1:b,3,1)



figure
poss=[min(LL_obs_stdnorm)-2:0.1:max(LL_obs_stdnorm)+2];

clear pdfmat pdfdat

for clust=1:selected_num_clusters
    pdfmat(:,clust)=MIXES(clust).*normpdf(poss,MEANS(clust),STDS(clust));
    pdfdat(:,clust)=MIXES(clust).*normpdf(LL_obs_stdnorm,MEANS(clust),STDS(clust));
    
end

plot(poss,pdfmat(:,1),'g')
hold on
if num_clusters>=2
    plot(poss,pdfmat(:,2),'r')
end

if num_clusters>=3
    plot(poss,pdfmat(:,3),'m')
    
    
end
plot(poss,sum(pdfmat,2),'k')
ksdensity(LL_obs_stdnorm)


[aa,idxx]=max(pdfdat');

figure
plot(LL_obs_stdnorm)
hold on
II=find(idxx==1);
plot(II,LL_obs_stdnorm(II),'ro','linewidth',[3]);

II=find(idxx==2);
plot(II,LL_obs_stdnorm(II),'go','linewidth',[3]);

II=find(idxx==3);
plot(II,LL_obs_stdnorm(II),'bo','linewidth',[3]);

II=find(idxx==4);
plot(II,LL_obs_stdnorm(II),'ko','linewidth',[3]);

II=find(idxx==5);
plot(II,LL_obs_stdnorm(II),'mo','linewidth',[3]);

if sum(idxx==2)>0
    legend('data','cluster 1','cluster 2')
end

if sum(idxx==2)>0&sum(idxx==3)>0
    legend('data','cluster 1','cluster 2','cluster 3')
end




