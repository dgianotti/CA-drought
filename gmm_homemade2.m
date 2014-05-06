clc;
clear;

load('CA_ids');

for i = 7:13
clearvars -except good_CA_IDs i
id = good_CA_IDs{i};


load(['LL_',id,'.mat']);

%load LL_049866


LLsim=LL_sim;
LLobs=LL_obs;
clear LL_obs LL_sim


LLsim=sum(LLsim,2)';
LLobs=sum(LLobs,2)';


SimLogLikelihood=-LLsim;
ObsLogLikelihood=-LLobs;

II=find(~isnan(ObsLogLikelihood));
ObsLogLikelihood=ObsLogLikelihood(II);




X=ObsLogLikelihood;


simtrace=SimLogLikelihood;


F=ksdensity(simtrace,X,'function','cdf');

XX=norminv(F,0,1);
%
% for ttt=1:500
%     ttt
% % test
% XI=simtrace(ceil(rand(length(X),1)*length(simtrace)));
%
% F=ksdensity(simtrace,XI,'function','cdf');
%
% XXI(:,ttt)=norminv(F,0,1);
% end

% test null


XX=norminv(F,0,1);





MU=0;
STD=1;
FORCEMEANSTD=1;



N=length(XX);

reps=1000; % max reps
LCRIT=0.001;




for clust=1:5
    clust
    
    
    
    for ens=1:100
        %1 randomly partition data into clusters
        idx=zeros(N,reps);
        idx(:,1)=ceil(rand(N,1)*clust);
        
        % estimate parameters assuming ident correct
        plotit=0;
        
        done=0;
        clear llsave
        replast = reps;
        
        for rep=1:reps
            
            if done==0
                        
                for ii=1:clust
                    if rep==1
                        II=find(idx(:,rep)==ii);
                        mixes(ii,rep)=length(II)./N;
                    end             
                    
                    if rep>1
            
                        pr=pdfdat(:,ii)./sum(pdfdat,2);
                        
                        mixes(ii,rep)=mean(pr);
                        
                        means(ii,rep)=sum(pr.*XX')./sum(pr);
                        
                        
%                         if ii==2
%                             means(ii,rep)=max(means(ii,rep),1); % at least one std to right
%                         end
%                         
%                         if ii==3
%                             means(ii,rep)=min(means(ii,rep),-1); % at least one std to left
%                         end
%                         
                        
                        stds(ii,rep)=sqrt(sum(pr.*(XX'-means(ii,rep)).^2)./sum(pr));
                        
                        stds(ii,rep)=max(stds(ii,rep),0.1);
                        
                        
                        if FORCEMEANSTD==1
                            if ii==1
                                means(ii,rep)=MU;
                                stds(ii,rep)=STD;
                                
                            end
                            
                        end
                        
                        
                    end
                    
                    % hard clust if rep==1
                    if rep==1
                        means(ii,rep)=mean(XX(II));
                        stds(ii,rep)=std(XX(II));
                    end
                    
                    
                end
                
                
                
                
                % create the prob model, and re-assign ident
                
                pdfdat=zeros(N,clust);
                for ii=1:clust
                    pdfdat(:,ii)=mixes(ii,rep).*normpdf(XX,means(ii,rep),stds(ii,rep));
                end
                
                [aa,idxx]=max(pdfdat');
                
                idx(:,rep+1)=idxx;
                
                if plotit==1
                    poss=[min(XX)-2:0.1:max(XX)+2];
                    
                    for ii=1:clust
                        pdfmat(:,ii)=mixes(ii,rep).*normpdf(poss,means(ii,rep),stds(ii,rep));
                    end
                    
                    plot(poss,pdfmat(:,1),'g')
                    hold on
                    plot(poss,pdfmat(:,2),'r')
                    
                    plot(poss,sum(pdfmat,2),'k')
                    ksdensity(XX)
                    
                end
                
                llsave(rep,1)=sum(log(sum(pdfdat,2)));
                
                if rep>1
                    if abs(llsave(rep,1)-llsave(rep-1,1))<LCRIT
                        done=1;
                        replast=rep;
                    end
                end
                
                
            end
            
        end
        
        
        aicsave(clust,ens)=-2*llsave(replast,1)+2*clust*2+2*(clust-1);
        loglikesave(clust,ens)=-2*llsave(replast,1);
        repsave(clust,ens)=replast;
        
        parsave(clust,1:clust,1,ens)=means(:,replast)';
        parsave(clust,1:clust,2,ens)=stds(:,replast)';
        parsave(clust,1:clust,3,ens)=mixes(:,replast)';
        
        
    end
end



for clust=1:5
    [a,b]=min(aicsave(clust,:));
    aicsave(clust,1)=a;
    loglikesave(clust,1)=loglikesave(clust,b);
    parsave(clust,1:clust,1,1)=parsave(clust,1:clust,1,b);
    parsave(clust,1:clust,2,1)=parsave(clust,1:clust,2,b);
    parsave(clust,1:clust,3,1)=parsave(clust,1:clust,3,b);
    
    
end



[a,b]=min(aicsave(:,1));

MEANS=parsave(b,1:b,1,1)
STDS=parsave(b,1:b,2,1)
MIXES=parsave(b,1:b,3,1)

clust=b;



%figure
XX(isinf(XX)) = NaN;
poss=[min(XX)-2:0.1:max(XX)+2];

clear pdfmat pdfdat

for ii=1:clust
    pdfmat(:,ii)=MIXES(ii).*normpdf(poss,MEANS(ii),STDS(ii));
    pdfdat(:,ii)=MIXES(ii).*normpdf(XX,MEANS(ii),STDS(ii));
    
end

plot(poss,pdfmat(:,1),'g')
hold on
if clust>=2
   plot(poss,pdfmat(:,2),'r')
end

if clust>=3
   plot(poss,pdfmat(:,3),'m')   
end
plot(poss,sum(pdfmat,2),'k')

print(gcf,'-dpng',['ClusterPDFs_',id,'.png']);

ksdensity(XX)










[aa,idxx]=max(pdfdat');


years = (2014-length(XX)):2013;

figure
plot(years,XX)
hold on
II=find(idxx==1);
h1 = plot(years(II),XX(II),'go','linewidth',3);

II=find(idxx==2);
h2 = plot(years(II),XX(II),'bo','linewidth',[3]);

II=find(idxx==3);
h3 = plot(years(II),XX(II),'ro','linewidth',[3]);

II=find(idxx==4);
h4 = plot(years(II),XX(II),'ko','linewidth',[3]);

II=find(idxx==5);
h5 = plot(years(II),XX(II),'mo','linewidth',[3]);


if sum(idxx==2)>0
    legend('data','Cluster 1','Cluster 2')
end

if sum(idxx==2)>0&sum(idxx==3)>0
    legend('data','Cluster 1','Cluster 2','Cluster 3')
end

ylim([-4,8]);

print(gcf,'-dpng',['YearlyClustering_',id,'.png']);
close all;
end

