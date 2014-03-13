function precip = get_precip_GSOD(datenums,lat,lon,max_distance)
% The function get_precip_GSOD takes as inputs:
% 
% datenums -- a vector of datenum values (see 'help datenum') corresponding
% to the days for which we need data,
%
% lat -- a single latitude value (+ for N, - for S)
% lon -- a single longitude value between -180 and 180 (- for W, + for E)
%
% max_distance -- a distance (in km) within which to search for data. If
% there is insufficient data within this range, NaNs will be returned for
% those days.

%datenums = missing_data_datenums;
%max_distance = 15;

% Make sure you have the GSOD station list:
ish_file_downloaded = (exist('ish-history.csv') == 2);
if ~ish_file_downloaded
    % Download:
    url = 'ftp://ftp.ncdc.noaa.gov/pub/data/gsod/ish-history.csv';
    urlwrite(url,'ish-history.csv');
end

fid = fopen('ish-history.csv');
fmt = '%s%s%s%*s%*s%*s%*s%s%s%*s%s%s';
F = textscan(fid, fmt, 'Delimiter',',', 'HeaderLines',1);
fclose(fid);

% F{1,1}: USAF ID
% F{1,2}: WBAN ID
% F{1,3}: STN Name
% F{1,4}: Lat
% F{1,5}: Lon
% F{1,6}: Start date
% F{1,7}: End date

% Remove all of the double quotes:
USAF = strrep(F{1,1}(:),'"','');
WBAN = strrep(F{1,2}(:),'"','');
STN_NAME = strrep(F{1,3}(:),'"','');

LAT = strrep(F{1,4},'"','');
LAT(strcmp(LAT,'')) = repmat({'-99999'},[sum(strcmp(LAT,'')),1]);
LAT = cellfun(@str2num,LAT);
LAT( ismember(LAT, [-9999,-99999,-999999]) ) = NaN;
LAT = LAT/1000;
LAT(LAT < -90) = NaN;
LAT(LAT > 90) = NaN;

LON = strrep(F{1,5},'"','');
LON(strcmp(LON,'')) = repmat({'-99999'},[sum(strcmp(LON,'')),1]);
LON = cellfun(@str2num,LON);
LON( ismember(LON, [-9999,-99999,-999999]) ) = NaN;
LON = LON/1000;
LON(LON < -180) = NaN;
LON(LON > 180) = NaN;

STARTDATE = strrep(F{1,6},'"','');
ENDDATE = strrep(F{1,7},'"','');

no_LAT = isnan(LAT);
no_LON = isnan(LON);
no_STARTDATE = strcmp('',STARTDATE);
%no_ENDDATE = strcmp('',ENDDATE);

bad_location = (no_LAT | no_LON | no_STARTDATE);


USAF(bad_location) = [];
WBAN(bad_location) = [];
STN_NAME(bad_location) = [];
LAT(bad_location) = [];
LON(bad_location) = [];
STARTDATE(bad_location) = [];
ENDDATE(bad_location) = [];

STARTYEAR = str2double(cellfun(@(x) x(1:4),STARTDATE,'un',0));
ENDYEAR = str2double(cellfun(@(x) x(1:4),ENDDATE,'un',0));

% Throw out data more than max_distance away:
distances = lldistkm([lat,lon],[LAT,LON]);
close_ones = distances <= max_distance;

USAF = USAF(close_ones);
WBAN = WBAN(close_ones);
STN_NAME = STN_NAME(close_ones);
LAT = LAT(close_ones);
LON = LON(close_ones);
STARTYEAR = STARTYEAR(close_ones);
ENDYEAR = ENDYEAR(close_ones);
distances = distances(close_ones);

% Sort them by proximity (closest first):
[distances,idx] = sort(distances);

USAF = USAF(idx);
WBAN = WBAN(idx);
STN_NAME = STN_NAME(idx);
LAT = LAT(idx);
LON = LON(idx);
STARTYEAR = STARTYEAR(idx);
ENDYEAR = ENDYEAR(idx);

precip = nan(size(datenums));

% For each station, download the data spanning datenums, and fill in any
% data you can:

for i = 1:length(distances)
    need_data_dates = datenums(isnan(precip));
    fprintf('%i days of data needed...\n',length(need_data_dates));
    fprintf('Trying station %i of %i nearby stations...\n',i,length(distances));
    
    min_year = str2double(datestr(min(need_data_dates),'yyyy')); % earliest year with no data
    max_year = str2double(datestr(max(need_data_dates),'yyyy')); % earliest year with no data

    % Download the data!
    % One year at a time
    for yr = min_year:max_year
        op_filename = [USAF{i},'-',WBAN{i},'-',num2str(yr),'.op'];
        gz_filename = [op_filename,'.gz'];
        
        url = ['ftp://ftp.ncdc.noaa.gov/pub/data/gsod/',...
            num2str(yr),'/', gz_filename];
        try
            if (~exist(['GSOD/',op_filename],'file') && (yr>=STARTYEAR{i}) && (yr <= ENDYEAR{i}))
                urlwrite(url,['GSOD/',gz_filename],'Timeout',20);
                gunzip(['GSOD/',gz_filename]);
                delete(['GSOD/',gz_filename]);
            end
        catch
            fprintf('No data found for station %s - %s for year %i...\n',...
                USAF{i},WBAN{i},yr);
        end
        
        if exist(['GSOD/',op_filename],'file') % if the file exists...
            fprintf('Filling year %i...\n',yr);
            
            fid = fopen(['GSOD/',op_filename]);
            % blah year, month, day, blah, precip, precip_flag
            fmt = ['%*14s','%4d','%2d','%2d','%*96s','%5s','%1s','%*14s'] ; 
            D = textscan(fid,fmt,'Whitespace','','Headerlines',1);
            fclose(fid);
            gsod_year = D{1};
            gsod_month = D{2};
            gsod_day = D{3};
            gsod_precip = cell2mat(cellfun(@str2num, D{4}, 'UniformOutput', false));
            gsod_precip( (gsod_precip == 99.99) & strcmp('I',D{5}) ) = nan; 

            % Now, assign any missing values with any data we just got!
            gsod_datenums = datenum( [num2str(gsod_year),...
                repmat('-',size(gsod_year)), num2str(gsod_month),...
                repmat('-',size(gsod_year)), num2str(gsod_day)] );
            
            % Loop over the missing data and fill if possible:
            for j = 1:length(need_data_dates)
                dnm = need_data_dates(j);
                if (ismember(dnm,gsod_datenums))
                    precip(datenums==dnm) = gsod_precip(gsod_datenums==dnm);
                end
            end
            
            if (numel(isnan(precip)) == 0)
                % Then we filled them all!
                fprintf('All precip data filled!\n');
                return;
            end
        end % if data for that year
           
    end % loop over years
    
end % loop over potential stations

if any(isnan(precip))
    fprintf(['We were not able to fill all of the missing data.\n',...
        'There are still %i days of missing data.\n'],sum(isnan(precip)));
end


end
