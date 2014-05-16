clear;
clc;



% Load station data:
for i = 1:13
    clearvars -except i % Because of all of the non-pre-allocated variables.  GAH!!
    
    maxclust=3;

    realrun=1;


addpath('C:/Users/gianotti/Documents/IntensityLib/');
load('CA_ids.mat');

id = good_CA_IDs{i};
load(sprintf('LL_%s.mat',id));

% Throw out 2014:
precip_obs(end,:) = [];

% Not sure about these...
Pobs = 365*nanmean(precip_obs,2)'; % Think this is right
POCCobs = nansum(precip_obs>0,2)'; % ?? THIS SHOULD PROBABLY BE NANMEAN TIMES THE NUMBER OF ~NANs, but not for testing in case it has to be an integer...


SimStn = load_stn_data(id,'SimStn');
Psim = sum(SimStn.intensity_data,2); % ??
POCCsim = sum(SimStn.intensity_data>0,2); % ??
N = length(Pobs); % number of obs years
n_sims = floor(length(Psim)/N);
Psim( (n_sims*N+1):end ) = []; % throw out extra sim years
Psim = Psim(:)';
POCCsim( (n_sims*N+1):end ) = []; % throw out extra sim years
POCCsim = POCCsim(:)';




if realrun==1
%    load test2
    
%    Psim=squeeze(sum(Pmat));
%    POCCsim=squeeze(sum(POCCmat));
%    Pobs=squeeze(sum(Pmatobs));
%    POCCobs=squeeze(sum(POCCmatobs));
 
    
    
%    Psim=reshape(Psim,1,length(Psim(:,1))*length(Psim(1,:)));
%    POCCsim=reshape(POCCsim,1,length(POCCsim(:,1))*length(POCCsim(1,:)));




X=Pobs;

simtrace=Psim;
F1=ksdensity(simtrace,X,'function','cdf');
Fsim1=ksdensity(simtrace,simtrace,'function','cdf');


X=POCCobs;

simtrace=POCCsim;
F2=ksdensity(simtrace,X,'function','cdf');
Fsim2=ksdensity(simtrace,simtrace,'function','cdf');


XX1=norminv(F1,0,1);
XX2=norminv(F2,0,1);

XX1sim=norminv(Fsim1,0,1);
XX2sim=norminv(Fsim2,0,1);


MU=0;
STD=1;
FORCEMEANSTD=1;

COVBASE=cov(XX1sim,XX2sim);



end




N=length(XX1);

reps=50; % max reps
LCRIT=0.001;

for clust=1:maxclust
    clust
    
    
    
    for ens=1:50
        %1 randomly partition data into clusters
        idx=zeros(N,reps);
        idx(:,1)=ceil(rand(N,1)*clust);
        
        % estimate parameters assuming ident correct
        plotit=0;
        
        done=0;
        clear llsave
           
          
            
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
                    
                    means1(ii,rep)=sum(pr.*XX1')./sum(pr);
                    means2(ii,rep)=sum(pr.*XX2')./sum(pr);
                   
                    
                
                    
                    covmat(1,1,ii,rep)=(sum(pr.*(XX1'-means1(ii,rep)).^2)./sum(pr));
                    covmat(2,2,ii,rep)=(sum(pr.*(XX2'-means2(ii,rep)).^2)./sum(pr));
                    covmat(1,2,ii,rep)=(sum(pr.*(XX1'-means1(ii,rep)).*(XX2'-means2(ii,rep)))./sum(pr));
                    covmat(2,1,ii,rep)=(sum(pr.*(XX2'-means2(ii,rep)).*(XX1'-means1(ii,rep)))./sum(pr));
                 
 %                   stds(ii,rep)=max(stds(ii,rep),0.1);
                    
                    
                    if FORCEMEANSTD==1
                        if ii==1
                            means(ii,rep)=MU;
                            stds(ii,rep)=STD;
                            covmat(:,:,ii,rep)=COVBASE;
                            
                        end
                        
                    end
                    
                    
                end
                
                % hard clust if rep==1
                if rep==1
                    means1(ii,rep)=mean(XX1(II));
                    means2(ii,rep)=mean(XX2(II));

                    yy=cov(XX1(II),XX2(II));
                    
                    covmat(:,:,ii,rep)=yy;
                end
                
                
            end
            
            
            
            
            % create the prob model, and re-assign ident
            
            pdfdat=zeros(N,clust);
            for ii=1:clust
                SIGMA=covmat(:,:,ii,rep);
                [~,p]=chol(SIGMA);
                if p>0
                                        
                    SIGMA(1,1)=SIGMA(1,1)+0.01;

                    SIGMA(2,2)=SIGMA(2,2)+0.01;
                end
                covmat(:,:,ii,rep)=SIGMA;
                pdfdat(:,ii)=mixes(ii,rep).*mvnpdf([XX1',XX2'],[means1(ii,rep),means2(ii,rep)],SIGMA);
            end
            
            [aa,idxx]=max(pdfdat');
            
            idx(:,rep+1)=idxx;
            

            llsave(rep,1)=sum(log(sum(pdfdat,2)));
            
            if rep>1
            if abs(llsave(rep,1)-llsave(rep-1,1))<LCRIT
                done=1;
                replast=rep;
            end
            end
            
            
            end
        
        end
        
        
       % aicsave(clust,ens)=-2*llsave(replast,1)+2*clust*2+2*(clust-1);
      %  for bivariate, more penalty parameters (means, mixes, covariances)
                aicsave(clust,ens)=-2*llsave(replast,1)+2*clust*4+2*(clust-1);

        loglikesave(clust,ens)=llsave(replast,1);
        repsave(clust,ens)=replast;
        
        parsave(clust,1:clust,1,ens)=means1(:,replast)';
        parsave(clust,1:clust,2,ens)=means2(:,replast)';
        parsave(clust,1:clust,3,ens)=mixes(:,replast)';
        
        for kk=1:clust
        parsave(clust,kk,4,ens)=covmat(1,1,kk,replast);
        parsave(clust,kk,5,ens)=covmat(1,2,kk,replast);
        parsave(clust,kk,6,ens)=covmat(2,1,kk,replast);
        parsave(clust,kk,7,ens)=covmat(2,2,kk,replast);
        end

        
    end
end

for clust=1:maxclust
    [a,b]=min(aicsave(clust,:));
    aicsave(clust,1)=a;
    loglikesave(clust,1)=loglikesave(clust,b);
    parsave(clust,1:clust,1,1)=parsave(clust,1:clust,1,b);
    parsave(clust,1:clust,2,1)=parsave(clust,1:clust,2,b);
    parsave(clust,1:clust,3,1)=parsave(clust,1:clust,3,b);
    parsave(clust,1:clust,4,1)=parsave(clust,1:clust,4,b);
    parsave(clust,1:clust,5,1)=parsave(clust,1:clust,5,b);
    parsave(clust,1:clust,6,1)=parsave(clust,1:clust,6,b);
    parsave(clust,1:clust,7,1)=parsave(clust,1:clust,7,b);
    
    
end


clear MEANS SIGMAS MIXES SimStn
[a,b]=min(aicsave(:,1));

for kk=1:b
MEANS(kk,:)=parsave(b,kk,1:2,1);
end


for kk=1:b
SIGMAS(1,1,kk)=(squeeze(parsave(b,kk,4,1)));
SIGMAS(1,2,kk)=(squeeze(parsave(b,kk,5,1)));
SIGMAS(2,1,kk)=(squeeze(parsave(b,kk,6,1)));
SIGMAS(2,2,kk)=(squeeze(parsave(b,kk,7,1)));

end


MIXES=parsave(b,1:b,3,1)

clust=b;




options = statset('Display','final');
% mu is k by d where k is mixes and d is dimenson;
% p is 1 by k
% sigma is d by d by k

obj=gmdistribution(MEANS,SIGMAS,MIXES);





%figure
%poss=[min(XX)-2:0.1:max(XX)+2];

clear pdfmat pdfdat

for ii=1:clust
    %pdfmat(:,ii)=MIXES(ii).*normpdf(poss,MEANS(ii),STDS(ii));
    %pdfdat(:,ii)=MIXES(ii).*normpdf(XX,MEANS(ii),STDS(ii));
                    
    pdfdat(:,ii)=MIXES(ii).*mvnpdf([XX1',XX2'],MEANS(ii,:),SIGMAS(:,:,ii));

end






[aa,idxx]=max(pdfdat');



figure
subplot(2,6,[1,2,7,8]);
h = ezcontour(@(x,y)pdf(obj,[x y]),[-8 6],[-8 6]);
hold on
for ii=1:clust
    II=find(idxx==ii);
if ii==1
plot(XX1(II),XX2(II),'ro','linewidth',[3]);
end
if ii==2
plot(XX1(II),XX2(II),'go','linewidth',[3]);
end
if ii==3
plot(XX1(II),XX2(II),'bo','linewidth',[3]);
end
end
title(id);


subplot(2,6,3:6);
plot(XX1)
hold on
II=find(idxx==1);
plot(II,XX1(II),'ro','linewidth',[3]);

II=find(idxx==2);
plot(II,XX1(II),'go','linewidth',[3]);

II=find(idxx==3);
plot(II,XX1(II),'bo','linewidth',[3]);

II=find(idxx==4);
plot(II,XX1(II),'ko','linewidth',[3]);

II=find(idxx==5);
plot(II,XX1(II),'mo','linewidth',[3]);

if sum(idxx==2)>0
    legend('data','cluster 1','cluster 2')
end

if sum(idxx==2)>0&sum(idxx==3)>0
    legend('data','cluster 1','cluster 2','cluster 3')
end





subplot(2,6,9:12);
plot(XX2)
hold on
II=find(idxx==1);
plot(II,XX2(II),'ro','linewidth',[3]);

II=find(idxx==2);
plot(II,XX2(II),'go','linewidth',[3]);

II=find(idxx==3);
plot(II,XX2(II),'bo','linewidth',[3]);

II=find(idxx==4);
plot(II,XX2(II),'ko','linewidth',[3]);

II=find(idxx==5);
plot(II,XX2(II),'mo','linewidth',[3]);

if sum(idxx==2)>0
    legend('data','cluster 1','cluster 2')
end

if sum(idxx==2)>0&sum(idxx==3)>0
    legend('data','cluster 1','cluster 2','cluster 3')
end
print(gcf,'-dpng',['plots/2d_',id,'.png']);



end % for loop over stns
