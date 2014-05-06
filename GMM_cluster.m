
% We want this to eventually be a function that takes LL_obs, LL_sim and
% returns the optimal GMM clustering

clear;
clc;
load('LL_046730.mat');

N = length(LL_obs);


% Sort the obs values:
obs_sorted = sort(LL_obs);


% Determine all possible GMM groupings with 1, 2, or 3 clusters (requiring
% monotonicity):
mu_hat = zeros(3,1);
sigma_hat = zeros(3,1);

clusters_best = ones(N,1);
k = 2; % just a single mean and variance

mu_hat(1) = mean(obs_sorted(clusters_best == 1));
sigma_hat(1) = std(obs_sorted(clusters_best == 1));
weights = [sum(clusters_best == 1), sum(clusters_best == 2), ...
    sum(clusters_best == 3)]/N;

L = nansum(log( weights(1)*normpdf(obs_sorted, mu_hat(1), sigma_hat(1))...
    + weights(2)*normpdf(obs_sorted, mu_hat(2), sigma_hat(2))...
    + weights(3)*normpdf(obs_sorted, mu_hat(3), sigma_hat(3)) ));


% BIC_best = -2*L + k*(log(N)+log(2*pi));
BIC_best = Inf;


% Determine the mean and variance of the sim cluster:
% For ALL possible clusters:
mu_hat(2) = mean(LL_sim(:));
sigma_hat(2) = std(LL_sim(:));

for num_1s = 0:(N-1)
    
    for num_2s = 2:(N-num_1s)
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
        
        L = nansum(log( weights(1)*normpdf(obs_sorted, mu_hat(1), sigma_hat(1))...
            + weights(2)*normpdf(obs_sorted, mu_hat(2), sigma_hat(2))...
            + weights(3)*normpdf(obs_sorted, mu_hat(3), sigma_hat(3)) ));
        
        BIC = -2*L + k*(log(N)+log(2*pi));
        
        if any([sum(clusters==1), sum(clusters==2), sum(clusters==3)]==1)
            % Require at least two in each cluster, problems otherwise
            BIC = Inf;
        end
        
        if BIC < BIC_best
            % Then it becomes the new best
            BIC_best = BIC;
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


hist(clusters_best)
figure;
ksdensity(LL_obs)
figure
x = floor(min([LL_sim(:);LL_obs(:)])):.01:ceil(max([LL_sim(:);LL_obs(:)]));
y = weights(1)*normpdf(x, mu_hat(1), sigma_hat(1))...
            + weights(2)*normpdf(x, mu_hat(2), sigma_hat(2))...
            + weights(3)*normpdf(x, mu_hat(3), sigma_hat(3));
plot(x,y);

% Unsort the obs values and cluster vector, and return the unsroted
% clusters, mean vector, std vector, and weights.



