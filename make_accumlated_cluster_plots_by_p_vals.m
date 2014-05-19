function make_accumlated_cluster_plots_by_p_vals(DOY)

load('CA_ids.mat');

marker_size = 40;

stn_names = {'Chula Vista', 'Colfax', 'Davis', 'Lemon Cove', 'Livermore',...    
	'Newport Beach','Ojai','Pasadena','Paso Robles','Petaluma',...
	'Redlands','Santa Barbara','Yreka'};

addpath('C:\Users\gianotti\Documents\GitHub\CA-drought\cbrewer\cbrewer\');

% Set the coloring intervals in terms of sigma values:
color_intervals = [ -Inf, -4.42, -3.89, -3.28, -2.56, 2.56, 3.29, 3.89, 4.42, Inf  ];

n_colors = length(color_intervals)-1;

cmap = cbrewer('div','RdYlBu',n_colors,'cubic');

for i = 1:length(good_CA_IDs)
    %fprintf('Clustering station %i of %i...\n',i,length(good_CA_IDs));
    stn_id = good_CA_IDs{i};
    % Load the data:
    filename = sprintf('LL_%s_accum_DOY%i.mat',stn_id,DOY);    
    load(filename);

    filename = sprintf('precip_%s_accum_DOY%i.mat',stn_id,DOY);
    load(filename);
    
	figure;
	for j = 1:9 % 9 plots -- {1,2,3} years, and {-LL, Total, Occ}
		switch j
			case 1 % 1yr LL
				data_obs = -LL_obs_1yr;
				data_sim = -LL_sim_1yr;
				N = length(data_obs);
				years = (2014-N):2013;
				make_title = true;
				title_text = sprintf('%s--Stn #%s: 1-yearly, Start DOY = %i',stn_names{i},stn_id,DOY);
				ylabel_text = '-LL';
				nino = nino_djfma_1yr;

			case 2 % 1yr Total
				data_obs = precip_obs_1yr;
				data_sim = precip_sim_1yr;
				N = length(data_obs);
				years = (2014-N):2013;
				make_title = false;
				ylabel_text = 'Total';
				nino = nino_djfma_1yr;
				
			case 3 % 1 yr Occ
				data_obs = occ_obs_1yr;
				data_sim = occ_sim_1yr;
				N = length(data_obs);
				years = (2014-N):2013;
				make_title = false;
				ylabel_text = 'Occ';
				nino = nino_djfma_1yr;
				
			case 4 % 2 yr LL
				data_obs = -LL_obs_2yr;
				data_sim = -LL_sim_2yr;
				N = length(data_obs);
				years = (2014.5-2*N):2:2012.5;
				make_title = true;
				title_text = '2-yearly';
				ylabel_text = '-LL';
				nino = nino_djfma_2yr;
				
			case 5 % 2 yr Total
				data_obs = precip_obs_2yr;
				data_sim = precip_sim_2yr;
				N = length(data_obs);
				years = (2014.5-2*N):2:2012.5;
				make_title = false;
				ylabel_text = 'Total';
				nino = nino_djfma_2yr;
				
			case 6 % 2 yr Occ
				data_obs = occ_obs_2yr;
				data_sim = occ_sim_2yr;
				N = length(data_obs);
				years = (2014.5-2*N):2:2012.5;
				make_title = false;
				ylabel_text = 'Occ';
				nino = nino_djfma_2yr;
				
			case 7 % 3 yr LL
				data_obs = -LL_obs_3yr;
				data_sim = -LL_sim_3yr;
				N = length(data_obs);
				years = (2015-3*N):3:2012;
				make_title = true;
				title_text = '3-yearly';
				ylabel_text = '-LL';
				nino = nino_djfma_3yr;
				
			case 8 % 3 yr Total
				data_obs = precip_obs_3yr;
				data_sim = precip_sim_3yr;
				N = length(data_obs);
				years = (2015-3*N):3:2012;
				make_title = false;
				ylabel_text = 'Total';
				nino = nino_djfma_3yr;
				
			case 9 % 3 yr Occ
				data_obs = occ_obs_3yr;
				data_sim = occ_sim_3yr;
				N = length(data_obs);
				years = (2015-3*N):3:2012;
				make_title = false;
				ylabel_text = 'Occ';
				nino = nino_djfma_3yr;
				
		end % switch		

		data_obs = (data_obs - mean(data_sim(:))) / std(data_sim(:));
		data_sim = (data_sim - mean(data_sim(:))) / std(data_sim(:));
		
		if (DOY > 1) % using something later than Jan1 start date
			years = years+1; % so that it ends with 2014
		end

		% Cluster data:
		%[clusters, mu_hat, sigma_hat, weights] = GMM_cluster(data_obs, data_sim);
		clusters = zeros(size(years));
		
		% plot time series
		subplot(9,7, (1:4) + 7*(j-1) );
		plot(years, data_obs,'-k');
		hold on;
		for col = 1:n_colors
			selected = (data_obs >= color_intervals(col)) & ...
				(data_obs <= color_intervals(col+1)); 
			scatter(years(selected), data_obs(selected), marker_size, cmap(col,:), 'filled', ...
				'MarkerEdgeColor','k');
			clusters(selected) = col;
		end % loop over colors
		
		ylm = get(gca,'Ylim');
		xlim([1900,2015]);

		if make_title
		    title(title_text);
		end
		ylabel(ylabel_text);
		
		% plot pdfs
		subplot(9,7,5 + 7*(j-1));
		y = linspace(ylm(1),ylm(2),200);
		plot( normpdf(y,mean(data_sim(:)), std(data_sim(:)) ), y, '-k');
		hold on;
		plot(ksdensity(data_obs,y),y,'-r');
		ylim(ylm);

		% plot QQ plot
		subplot(9, 7, 6 + 7*(j-1));
		qqplot(data_sim(:),data_obs(:));
		xlabel('Sim Quantile');
		ylabel('Obs Quantile');
		ylim(ylm);
		
		% Plot ENSO correlation
		subplot(9, 7, 7 + 7*(j-1));
		for col = 1:n_colors
			clust = clusters == col;
			scatter(nino(clust), data_obs(clust), marker_size, cmap(col,:), 'filled', ...
				'MarkerEdgeColor','k');
			hold on;
		end % loop over colors
		nino_r2 = corrcoef(data_obs,nino);
		xlabel(sprintf('Nino 3.4 Dec-Apr: R^2 = %1.2f',nino_r2(2)^2));
		

		
	end % Loop over plot rows
	
	set(gcf,'Position',[1 1 1680 1800],...
		'PaperSize',[15,18], 'PaperPosition',[1,1,14,17]);
	outfilename = sprintf('plots/Pval_TS_%s_DOY%i.pdf',stn_id,DOY);
    print(gcf,'-dpdf',outfilename);
 
end % loop over stations
    
close all;
end % function