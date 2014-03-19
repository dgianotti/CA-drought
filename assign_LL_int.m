function LL_int = assign_LL_int(int, params)
% The function assign_LL_int calculates the log likelihood (really log
% probability) of getting the intensity values given in the int vector,
% given the parameters in params.
%
% int is an N x 1 vector of precipitation intensities. 
%
% params is a 5 x 1 vector of parameters for our sum of two gamma
% distributions model. The model looks like:
% [k1; theta1; k2; theta2; weight]
% with pdf = params(5)*gampdf(data,params(1),params(2)) ...
%              + (1-params(5))*gampdf(data,params(3),params(4))

% We're going to define a probability of witnessing a certain intensity by
% defining bins and seeing the probability of falling into those bins.
% There are three ways I can think of to define bin sizes:
% 1) as simple constant increments (like 0-1mm, 1-2mm...)
% 2) as a fraction of the distribution parameters (like mean/10, 2*mean/10,
% 3*mean/10...)
% 3) a fixed number, with the range defined somehow by the distribution
% (100 increments--evenly spaced in depth--between 0 and 4sigma above the mean).

% To deal with missing data, we will first assign all nans to be 1,
% determine the individual LL, then reassign nans at the end:
missing_data = isnan(int);
int(missing_data) = 1;

addpath('C:\Users\gianotti\Documents\IntensityLib');
% We'll try 0.254mm increments starting at 0:
increment_size = 0.254;

% maximum foreseeable daily precip? ~1.5m?
max_precip = floor(1500/increment_size)*increment_size;

bin_edges = 0:increment_size:max_precip;

cdfs = mixgamcdf(bin_edges(2:end),params); % m x 1
cdfs = [0,cdfs]; % just in case;

bin_probabilities = diff(cdfs); % these should be all of the probabilities.

% The sum of bin_probabilities should be 1. We'll make sure:
bin_probabilities(end) = 1-sum(bin_probabilities(1:(end-1))); % now sum to 1
% And then  we'll fix the end of the cdfs too:
cdfs(end) = 1;

% Now we need to figure out which bin each of the values in int fit into.
[~,which_bin] = histc(int,bin_edges);

p_int = bin_probabilities(which_bin);
p_int(int == 0) = 1; % so that LL = 0
LL_int = log(p_int)';

% Now re-assign the missing data to nans:
LL_int(missing_data) = nan;
end


