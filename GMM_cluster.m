function [clusters, mu_hat, sigma_hat, weights] = GMM_cluster(obs, sim)
% The function GMM_cluster takes two arguments, obs and sim, and
% clusters the obs values according to a Gaussian mixture model with
% between 1-3 components (selected by BIC). 

if ~isvector(obs)
    error('LL_obs should be a vector (hopefully with about 100 elements). Aborting!');
end

N = length(obs);

% Sort the obs values:
[obs_sorted, sort_idx] = sortrows(obs(:));

% Determine all possible GMM groupings with 1, 2, or 3 clusters (requiring
% monotonicity):
mu_hat = zeros(3,1);
sigma_hat = zeros(3,1);

clusters_best = ones(N,1);
%k = 2; % just a single mean and variance

%mu_hat(1) = mean(obs_sorted(clusters_best == 1));
%sigma_hat(1) = std(obs_sorted(clusters_best == 1));
%weights = [sum(clusters_best == 1), sum(clusters_best == 2), ...
%    sum(clusters_best == 3)]/N;

%L = nansum(log( weights(1)*normpdf(obs_sorted, mu_hat(1), sigma_hat(1))...
%    + weights(2)*normpdf(obs_sorted, mu_hat(2), sigma_hat(2))...
%    + weights(3)*normpdf(obs_sorted, mu_hat(3), sigma_hat(3)) ));


% BIC_best = -2*L + k*(log(N)+log(2*pi));
% BIC_best = Inf;
AIC_best = Inf;


% Determine the mean and variance of the sim cluster:
% For ALL possible clusters:
mu_hat(2) = mean(sim(:));
sigma_hat(2) = std(sim(:));

for num_1s = 0:(N-1)
    
    % Require at least N/4 values to be in the central cluster
    for num_2s = ceil(N/4):(N-num_1s)
        clusters = [ones(num_1s,1); repmat(2,[num_2s,1]); ...
            repmat(3,[N-num_1s-num_2s,1])];
        
        % Need a mean and a variance for each cluster, so
        k = 2*numel(unique(clusters));
        
        % Calculate estimates of mu, sigma:
        mu_hat(1) = mean(obs_sorted(clusters == 1));
        sigma_hat(1) = std(obs_sorted(clusters == 1));
        %mu_hat(2) = mean(obs_sorted(clusters == 2));
        %sigma_hat(2) = std(obs_sorted(clusters == 2));
        mu_hat(3) = mean(obs_sorted(clusters == 3));
        sigma_hat(3) = std(obs_sorted(clusters == 3));
        
        weights = [sum(clusters == 1), sum(clusters == 2), ...
            sum(clusters == 3)]/N;
        
        mu_hat(isnan(mu_hat)) = 0;
        sigma_hat(isnan(sigma_hat)) = 1; 
        
        % Set the minimum std for clusters 1 and 3 to be the std 
        % of the smallest/largest values so that the std is never 0
        % (leads to BIC problems):
        if (sum(clusters==1) == 1)
            sigma_hat(1) = max(sigma_hat(1),std(obs_sorted(1:3)));
        end
        if (sum(clusters==3) == 1)
            sigma_hat(3) = max(sigma_hat(3),std(obs_sorted((end-2):end)));
        end
        
        
        L = nansum(log( weights(1)*normpdf(obs_sorted, mu_hat(1), sigma_hat(1))...
            + weights(2)*normpdf(obs_sorted, mu_hat(2), sigma_hat(2))...
            + weights(3)*normpdf(obs_sorted, mu_hat(3), sigma_hat(3)) ));
        
        % BIC = -2*L + k*(log(N)+log(2*pi));
        AIC = -2*L + 2*k*N/(N-k-1);
        
        % I think this is solved by our minimum variance requirement:
%         if any([sum(clusters==1), sum(clusters==2), sum(clusters==3)]==1)
%             % Require at least two in each cluster, problems otherwise
%             BIC = Inf;
%         end
        
        if AIC < AIC_best
            % Then it becomes the new best
            AIC_best = AIC;
            clusters_best = clusters;
        end    
    end
end

mu_hat(1) = mean(obs_sorted(clusters_best == 1));
sigma_hat(1) = std(obs_sorted(clusters_best == 1));
mu_hat(3) = mean(obs_sorted(clusters_best == 3));
sigma_hat(3) = std(obs_sorted(clusters_best == 3));

mu_hat(isnan(mu_hat)) = 0;
sigma_hat(isnan(sigma_hat)) = 1;

weights = [sum(clusters_best == 1), sum(clusters_best == 2), ...
    sum(clusters_best == 3)]/N;


% figure;
% subplot(2,3,[1,4]);
% hist(clusters_best)
% title('Histogram of cluster assignments');
% subplot(2,3,2:3);
% ksdensity(LL_obs);
% title('ksdensity of LL_{obs}');
% subplot(2,3,5:6);
% x = floor(min([LL_sim(:);LL_obs(:)])):.01:ceil(max([LL_sim(:);LL_obs(:)]));
% y = weights(1)*normpdf(x, mu_hat(1), sigma_hat(1))...
%             + weights(2)*normpdf(x, mu_hat(2), sigma_hat(2))...
%             + weights(3)*normpdf(x, mu_hat(3), sigma_hat(3));
% plot(x,y);
% title('Assigned cluster PDFs');


% Unsort the obs values and cluster vector, and return the unsorted
% clusters, mean vector, std vector, and weights.

tmp = sortrows([sort_idx,clusters_best(:)],1);
clusters = tmp(:,2);

end % function

