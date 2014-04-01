function [std_norm_pdf_vals, uniform_pdf_vals, ranks] = empirical_2_normal_via_order_stats(X)
% The function empirical_2_normal_via_order_stats takes as input a vector
% or matrix X and returns the probability densities of those values as they
% are re-assigned to a standard normal distribution via the order
% statistics and an inverse cdf transform (a histogram of std_norm_pdf_vals
% will essential be perfectly gaussian, assigned as the most likely values
% of the order statistics).
%

[J,K] = size(X);
 
N = J*K; % total number of entries

X = X(:); % Work with it as a vector, then reshape at the end

[X_sorted,rank] = sort(X);

order_stats_sorted = (1:N)/(N+1); % These are the most likely values for a uniform variable
uniform_pdf_vals = reshape(order_stats_sorted(rank),[J,K]);

std_norm_pdf_vals = norminv(uniform_pdf_vals,0,1); 


rank = reshape(rank,[J,K]);
 