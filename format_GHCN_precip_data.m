function [datenums, precip] = format_GHCN_precip_data(filename,varargin)
% The function format_GHCN_precip_data(filename) takes as its input a single file
% of daily data from the GHCN/USHCN (usually .dly) and outputs two vectors:
%
% datenums is an N x 1 vector of date-numbers (see 'help datenum'),
% 
% precip is an N x 1 vector of precipitation amounts in mm corresponding to
% the date in datenums.
%
%
% format_GHCN_precip_data(filename,'PadFirstLastYears',true) extends
% datenums to Jan 1 of the first year in the dataset and to Dec 31 of the
% last year of the dataset  [default: false]. Corresponding precip values
% are NaN.
%
%
% format_GHCN_precip_data(filename,'ExcludeLeapDays',true) removes all
% February 29ths from the data set [default: false].

ExcludeLeapDays = false; % default
PadFirstLastYears = false; % default

if nargin == 3 % Flags for ExcludeLeapDays or PadFirstLastYears
    if (strcmpi(varargin{1},'ExcludeLeapDays'))
        ExcludeLeapDays = varargin{2};
    elseif (strcmpi(varargin{1},'PadFirstLastYears'))
        PadFirstLastYears = varargin{2};
    else
        error('The second argument to format_GHCN_precip_data should be either ExcludeLeapDays or PadFirstLastYears. Aborting!\n');
    end
elseif nargin == 5 % Flags for ExcludeLeapDays AND PadFirstLastYears
    if (strcmpi(varargin{1},'ExcludeLeapDays'))
        ExcludeLeapDays = varargin{2};
    elseif (strcmpi(varargin{1},'PadFirstLastYears'))
        PadFirstLastYears = varargin{2};
    else
        error('The second argument to format_GHCN_precip_data should be either ExcludeLeapDays or PadFirstLastYears. Aborting!\n');
    end
    
    if (strcmpi(varargin{3},'ExcludeLeapDays'))
        ExcludeLeapDays = varargin{4};
    elseif (strcmpi(varargin{3},'PadFirstLastYears'))
        PadFirstLastYears = varargin{4};
    else
        error('The fourth argument to format_GHCN_precip_data should be either ExcludeLeapDays or PadFirstLastYears. Aborting!\n');
    end
    
elseif nargin ~= 1 
    error('The function format_GHCN_precip_data should have 1, 3, or 5 inputs. Aborting!\n');
end
    

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

% We need different sized vectors if we're keeping leap days:
if ExcludeLeapDays
   year_vec = zeros(range(YEAR)*365,1);
else
    % Add in the necessary number of leap days:
    yrs = min(YEAR):max(YEAR);
    num_leap_days = sum(mod(yrs,4)==0);
    year_vec = zeros(range(YEAR)*365+num_leap_days,1);
end
month_vec = year_vec;
day_vec = year_vec;
precip = year_vec;

n_rows = size(YEAR,1);
idx = 1;

num_months = 12*range(YEAR);
for i = 1:num_months
    yr = min(YEAR) + ceil((i)/12) - 1;
    mnth = mod(i-1,12)+1;
    
    if (~ExcludeLeapDays && (mod(yr,4) == 0)) % If a leap year and we're keeping leap days
        days_per_month(2) = 29;
    else 
        days_per_month(2) = 28;
    end
    
    n_days = days_per_month(mnth);
    year_vec(idx:(idx+n_days-1)) = yr;
    month_vec(idx:(idx+n_days-1)) = mnth;
    day_vec(idx:(idx+n_days-1)) = 1:n_days;
    
    % Now find the right row:
    row_match = (YEAR == yr) & (MONTH == mnth);
    if (sum(row_match) > 1) % That's bad... multiple copies of the same month/year in the file?
        fprintf('There seem to be multiple copies of the precip data for year %i, month %i.\n',yr,mnth);
        error('Aborting!\n');
    elseif (sum(row_match) == 0)
        % No data for that month, enter NaNs:
        precip(idx:(idx+n_days-1)) = NaN;
    else % Data!
        precip(idx:(idx+n_days-1)) = VALUE(row_match,1:n_days);
    end    
    idx = idx + n_days;    
end

precip(precip < 0) = NaN;

precip = precip/10; % Now in mm

datenums = datenum( [num2str(year_vec),repmat('-',[length(precip),1]),...
    num2str(month_vec,'%02d'),repmat('-',[length(precip),1]),num2str(day_vec,'%02d')] );

if ~PadFirstLastYears
    % Remove the leading and trailing NaNs:
    first_data_idx = find(~isnan(precip),1);
    last_data_idx = find(~isnan(precip),1,'last');
    
    to_remove = [1:(first_data_idx-1),(last_data_idx+1):length(year_vec)];
    datenums(to_remove) = [];
    precip(to_remove) = [];
    
end


end