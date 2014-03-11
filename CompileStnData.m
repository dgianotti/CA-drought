% Let's get the (current) data for the CA stations!

clear;
clc;

load('CA_ids.mat');

% for i = 1:length(CA_IDs)
id = CA_IDs{1}; % Berkeley
A = importdata('301137.csv',',');
dates = A.data(:,1);
precip = A.data(:,2)/10; % in mm

% throw out the data before 2010:
old_data = dates<20100000;
dates(old_data) = [];
precip(old_data) = [];

plot(dates,precip)
% Well... that doesn't look good...