if 1==2
    
% load 082944_SON_ObsSimLL
% SimLogLikelihood(4,537) = -178.5;

m1=mean(mean(SimLogLikelihood));
sig1=mean(std(SimLogLikelihood'));
poss=[m1-5*sig1:1:m1+5*sig1];






% first try for sim data where should get one cluster only

options = statset('MaxIter',1000,'Display','Final');
replicates=10;
X=ObsLogLikelihood;
clear aicsave_sim bicsave_sim

for mmm=1:50
    mmm
X=SimLogLikelihood(:,mmm);

for nnn=1:3
%gm = gmdistribution.fit(X,nnn,'regularize',[0.01],'Options',options,'CovType','full','Replicates',replicates);
gm = gmdistribution.fit(X,nnn,'regularize',[0.01],'Options',options,'CovType','diag','Replicates',replicates);

aicsave_sim(nnn,mmm)=gm.AIC;
bicsave_sim(nnn,mmm)=gm.BIC;

end

end


[a,bchoiceaic]=min(aicsave_sim);

[a,bchoicebic]=min(bicsave_sim);
figure

hist([bchoiceaic',bchoicebic'])
title('number of times each order chosen when applied to null (simulated likelihoods)')
legend('AIC','BIC')



X=ObsLogLikelihood;

clear aicsave bicsave
replicates=50;

for nnn=1:3
%gm = gmdistribution.fit(X,nnn,'regularize',[0.01],'Options',options,'CovType','full','Replicates',replicates);
gm = gmdistribution.fit(X,nnn,'regularize',[0.01],'Options',options,'CovType','diag','Replicates',replicates);

aicsave(nnn,1)=gm.AIC;
bicsave(nnn,1)=gm.BIC;

end





[aa,bb]=min(bicsave);
nn=bb

figure
plot(bicsave)
hold on
plot(bb,bicsave(bb),'ko','linewidth',[3])
title('BIC for observed data showing optimal number of clusters')
ylabel('BIC')
xlabel('number of clusters')


% redo fit for chosen number of clusters

gm = gmdistribution.fit(X,bb,'Options',options,'CovType','full','Replicates',replicates);
gm = gmdistribution.fit(X,bb,'regularize',[0.01],'Options',options,'CovType','diag','Replicates',replicates);

gm.NlogL;

gm.PComponents
gm.mu
gm.Sigma
clear pdfmat
for comp=1:bb
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


end





% home made version

% load 082944_SON_ObsSimLL
% SimLogLikelihood(4,537) = -178.5;

load 353770_likelihood_data

m1=mean(mean(SimLogLikelihood));
sig1=mean(std(SimLogLikelihood'));
poss=[m1-5*sig1:1:m1+5*sig1];

X=ObsLogLikelihood;





mu1=median(X);
mu1=mean(X);

sig1=mean(std(SimLogLikelihood'));


lowendmu=(mu1-4*sig1):(3*sig1/10):(mu1-1*sig1);

highendmu=(mu1+1*sig1):(3*sig1/10):(mu1+4*sig1);

lowendsig=sig1.*[0.1:0.025:0.2];
highendsig=sig1.*[0.1:0.025:0.2];

mix1s=[0.8:0.025:1];
mix2s=[0:0.025:0.1];


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

[aa,bb]=sort(-2*llsave(1:cc)+2*pensave(1:cc));
    

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



plot(poss,pdf1,'g','linewidth',[3])
hold on
plot(poss,pdf2,'r','linewidth',[3])
plot(poss,pdf3,'b','linewidth',[3])
plot(poss,pdfall,'k','linewidth',[3])




% now for the associateion part...


pdf1dat=mix1.*normpdf(X,mu1,sig1);
pdf2dat=mix2.*normpdf(X,mu2,sig2);
pdf3dat=+mix3.*normpdf(X,mu3,sig3);

[aa,idx]=max([pdf1dat,pdf2dat,pdf3dat]');

% idx has in it the cluster of each data point



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

if sum(idx==2)>0
legend('data','cluster 1','cluster 2')
end

if sum(idx==2)>0&sum(idx==3)>0
legend('data','cluster 1','cluster 2','cluster 3')
end



