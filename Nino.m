clear;
clc;

load('nino34.mat');

load('CA_ids.mat');
for i = 1:13
    id = good_CA_IDs{i};
    filename = sprintf('LL_%s_accum_DOY1.mat',id);
    load(filename);
    filename = sprintf('precip_%s_accum_DOY%i.mat',id,1);
    load(filename);
    
    dataLL = LL_obs_1yr_stdnorm(1:(end-1));
    dataprecip = precip_obs_1yr_stdnorm(1:(end-1));
    dataocc = occ_obs_1yr_stdnorm(1:(end-1));
    
    N = length(dataLL);
    
    nino = nino34_annual( (end-N):(end-1));
 
    figure;
    subplot(3,1,1);
    plot(nino,-dataLL,'o','LineStyle','none')
    title(id);
    ylabel('-LL');
    
    subplot(3,1,2);
    plot(nino,dataprecip,'o','LineStyle','none')
    ylabel('Total Precip');
    
    subplot(3,1,3);
    plot(nino,dataocc,'o','LineStyle','none')
    ylabel('Occurrence');
    xlabel('Nino 3.4');
    
    print(gcf,'-dpng',['Nino34_',id,'.png']);
    
end
    