function [datenums, precip] = format_GHCN_precip_data(filename)
% The function format_GHCN_precip_data takes as its input a single file
% of daily data from the GHCN/USHCN (usually .dly) and outputs two vectors:
%
% datenums is an N x 1 vector of date-numbers (see 'help datenum')
fid = fopen(filename);
% ------------------------------
% Variable   Columns   Type
% ------------------------------
% ID            1-11   Character
% YEAR         12-15   Integer
% MONTH        16-17   Integer
% ELEMENT      18-21   Character
% VALUE1       22-26   Integer
% MFLAG1       27-27   Character
% QFLAG1       28-28   Character
% SFLAG1       29-29   Character
% VALUE2       30-34   Integer
% MFLAG2       35-35   Character
% QFLAG2       36-36   Character
% SFLAG2       37-37   Character
%   .           .          .
%   .           .          .
%   .           .          .
% VALUE31    262-266   Integer
% MFLAG31    267-267   Character
% QFLAG31    268-268   Character
% SFLAG31    269-269   Character
% ------------------------------
fmt = ['%11s %4s %2s %4s ',...
    '%5s %1s %1s %1s %5s %1s %1s %1s %5s %1s %1s %1s %5s %1s %1s %1s ',...
    '%5s %1s %1s %1s %5s %1s %1s %1s %5s %1s %1s %1s %5s %1s %1s %1s ',...
    '%5s %1s %1s %1s %5s %1s %1s %1s %5s %1s %1s %1s %5s %1s %1s %1s ',...
    '%5s %1s %1s %1s %5s %1s %1s %1s %5s %1s %1s %1s %5s %1s %1s %1s ',...
    '%5s %1s %1s %1s %5s %1s %1s %1s %5s %1s %1s %1s %5s %1s %1s %1s ',...
    '%5s %1s %1s %1s %5s %1s %1s %1s %5s %1s %1s %1s %5s %1s %1s %1s ',...
    '%5s %1s %1s %1s %5s %1s %1s %1s %5s %1s %1s %1s %5s %1s %1s %1s ',...
    '%5s %1s %1s %1s %5s %1s %1s %1s %5s %1s %1s %1s'];                    
D = textscan(fid,fmt,'Whitespace','');
fclose(fid);
% Now convert the numeric data to numeric:

D{2} = cellfun(@str2num, D{2}, 'UniformOutput', false);
YEAR = cell2mat(D{2});

D{3} = cellfun(@str2num, D{3}, 'UniformOutput', false);
MONTH = cell2mat(D{3});

VALUE = nan(size(YEAR),31);
day = 1;
for i = 5:4:125 % All of the 'VALUE' columns
     D{i} = cellfun(@str2num, D{i}, 'UniformOutput', false);
     VALUE(:,day) = cell2mat(D{i});
     day = day+1;
end

% Boy, now we have a cell array of 128 columns...
% Let's get just the precip data:
precip_rows = strcmp(D{4},'PRCP');

YEAR = YEAR(precip_rows,:);
MONTH = MONTH(precip_rows,:);
VALUE = VALUE(precip_rows,:);

days_per_month = [31,28,31,30,31,30,31,31,30,31,30,31]; % Ignoring leap days

% Now, for each row, we'll add the dates to the vector of dates and the
% precip to the vector of precip:
year_vec = zeros(range(YEAR)*365,1);
month_vec = year_vec;
day_vec = year_vec;
precip_vec = year_vec;

n_rows = size(YEAR,1);
idx = 1;
for i = 1:n_rows
    n_days = days_per_month(MONTH(i));
    year_vec(idx:(idx+n_days-1)) = repmat(YEAR(i),[n_days,1]);
    month_vec(idx:(idx+n_days-1)) = repmat(MONTH(i),[n_days,1]);
    day_vec(idx:(idx+n_days-1)) = 1:n_days;
    precip_vec(idx:(idx+n_days-1)) = VALUE(i,1:n_days);
    
    idx = idx + n_days;    
end

precip_vec(precip_vec < 0) = nan;

precip_vec = precip_vec/10; % Now in mm

unfilled = (month_vec == 0);
year_vec(unfilled) = [];
month_vec(unfilled) = [];
day_vec(unfilled) = [];
precip_vec(unfilled) = [];


date_vec = datenum( [num2str(year_vec),repmat('-',[length(precip_vec),1]),...
    num2str(month_vec,'%02d'),repmat('-',[length(precip_vec),1]),num2str(day_vec,'%02d')] );


end