% Let's get the (current) data for the CA stations!

clear;
clc;

load('CA_ids.mat');

%% Download data:
% for i = 1:length(CA_IDs)
id = CA_IDs{1}; % Berkeley
% Download the latest USHCN data:
url = ['http://www1.ncdc.noaa.gov/pub/data/ghcn/daily/hcn/USC00',id,'.dly'];
filename = [id,'.dly'];
urlwrite(url,filename);

%% Read in fixed-width data:

[datenums,precip] = format_GHCN_precip_data(filename,'PadFirstLastYears',true,'ExcludeLeapDays',true);


%%

A = importdata('301137.csv',',');
dates = A.data(:,1);
precip = A.data(:,2)/10; % in mm

% throw out the data before 2010:
old_data = dates<20100000;
dates(old_data) = [];
precip(old_data) = [];

%% Convert dates to date numbers:
% Convert dates to strings:
dates = num2str(dates);
n = size(dates,1);

% First 4 digits are yyyy
yyyy = dates(:,1:4);
mm = dates(:,5:6);
dd = dates(:,7:8);
date_str = horzcat(yyyy,repmat('-',[n,1]),mm,repmat('-',[n,1]),dd); % now yyyy-mm-dd
date_nums = datenum(date_str);

%%
%plot(precip)
plot(date_nums, precip)
datetick('x');
% Well... that doesn't look good...