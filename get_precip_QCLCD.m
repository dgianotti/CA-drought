function precip = get_precip_QCLCD(datenums,lat,lon)
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

precip = nan(size(datenums));

fprintf('%i days of data needed...\n',length(datenums));

min_year = str2double(datestr(min(datenums),'yyyy')); % earliest year with no data
max_year = str2double(datestr(max(datenums),'yyyy')); % latest year with no data

datenums_yyyy = cellstr(datestr(datenums,'yyyy'));
datenums_mm = cellstr(datestr(datenums,'mm'));
datenums_dd = cellstr(datestr(datenums,'dd'));


for yr = min_year:max_year
    
    for mnth = 1:12
        YYYYMM = [num2str(yr),num2str(mnth,'%02i')];
        % If there are any missing dates in this month, get the data and
        % fill in precip:
        if any(strcmp(datenums_yyyy,num2str(yr)) &...
                strcmp(datenums_mm,num2str(mnth,'%02i')))
            
            % Check to see if data is downloaded:
            if ~exist(['QCLCD/',YYYYMM,'station.txt'],'file') % if the station file doesn't exist
                % Download it:
                filename = ['QCLCD',YYYYMM,'.zip'];
                url = ['http://cdo.ncdc.noaa.gov/qclcd_ascii/',filename];
                fprintf('Downloading file %s ... \n',filename);
                urlwrite(url,['QCLCD/',filename]);
                fprintf('Unzipping file %s ... \n',filename);
                unzip(['QCLCD/',filename],'QCLCD');
                delete(['QCLCD/',filename]); 
                delete(['QCLCD/',YYYYMM,'hourly.txt']); 
                delete(['QCLCD/',YYYYMM,'monthly.txt']); 
                delete(['QCLCD/',YYYYMM,'precip.txt']); 
                delete(['QCLCD/',YYYYMM,'remarks.txt']);                 
            end
  
            % Find the closest station:
            fid = fopen(['QCLCD/',YYYYMM,'station.txt']);
            fmt = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s';
            F = textscan(fid, fmt, 'Delimiter','|', 'HeaderLines',1);
            fclose(fid);
            
            % F{1}: WBAN
            % F{2}: WMO
            % F{3}: CallSign
            % F{4}: ClimateDivisionCode
            % F{5}: ClimateDivisionStateCode
            % F{6}: ClimateDivisionStationCode
            % F{7}: Name
            % F{8}: State
            % F{9}: Location
            % F{10}: Latitude
            % F{11}: Longitude
            % F{12}: GroundHeight
            % F{13}: StationHeight
            % F{14}: Barometer
            % F{15}: TimeZone
            
            % We just need the ID, lat, and lon:
            WBAN = F{1};
            LAT = cellfun(@str2num,F{10});
            LON = cellfun(@str2num,F{11});
            
            % Find the closest station:
            distances = lldistkm([lat,lon],[LAT,LON]);
            [d,idx] = min(distances);
            wban_id = WBAN{idx};

            seemingly_bad_stations = {'23239'};
            while any(ismember(seemingly_bad_stations, wban_id))
                distances(idx) = inf;
                [d,idx] = min(distances);
                wban_id = WBAN{idx};
            end
                
            loc = cell2mat(F{9}(idx));
            fprintf('The closest station is %1.1f km away at %s,\nWBAN ID: %s.\n',...
                d,loc,wban_id);
            
                

            % Now fill in the missing data with data from that station:
            clear('F');
            fid = fopen(['QCLCD/',YYYYMM,'daily.txt']);
            fmt = ['%s%s%s%s%s%s%s%s%s%s',...
                '%s%s%s%s%s%s%s%s%s%s',...
                '%s%s%s%s%s%s%s%s%s%s',...
                '%s%s%s%s%s%s%s%s%s%s',...
                '%s%s%s%s%s%s%s%s%s%s'];
        
            F = textscan(fid, fmt, 'Delimiter',',', 'HeaderLines',1);
            fclose(fid);
            
                       
            QCLCD_precip = str2double(F{31});
            % Set "trace" amounts to zero:
            QCLCD_precip(isnan(QCLCD_precip)) = 0;
            
            QCLCD_date = F{2};
            QCLCD_wban = F{1};
            
            % For each needed data point in that month, fill in with the
            % proper precip value:
            this_months_datenums_dd = cell2mat(datenums_dd(strcmp(datenums_yyyy,num2str(yr)) &...
                strcmp(datenums_mm,num2str(mnth,'%02i'))));
            this_months_date_strings = [repmat(YYYYMM,...
                [size(this_months_datenums_dd,1),1]),this_months_datenums_dd];

            for i = 1:size(this_months_date_strings,1)
                % For each missing date:
                YYYYMMDD = this_months_date_strings(i,:);
                if ( sum(strcmp(YYYYMMDD,QCLCD_date) & ...
                        strcmp(wban_id,QCLCD_wban)) == 1)
                
                    precip( datenums == datenum(YYYYMMDD,'yyyymmdd') ) = ...
                        QCLCD_precip( strcmp(YYYYMMDD,QCLCD_date) & ...
                        strcmp(wban_id,QCLCD_wban) );
                else
                    fprintf('There are %i values for the date %s.\n Cannot fill data!\n',...
                        sum(strcmp(YYYYMMDD,QCLCD_date) & ...
                        strcmp(wban_id,QCLCD_wban)),YYYYMMDD);
                end
            end
            
            
            if (sum(isnan(precip)) == 0)
                % Then we filled them all!
                fprintf('All precip data filled!\n');
                return;
            end
            
        end % If statement for missing data this month
        fprintf('Filled %i days of data, %i days remaining...\n',...
            size(this_months_datenums_dd,1),sum(isnan(precip)));
        
    end % loop over months
end % loop over years

if any(isnan(precip))
    fprintf(['We were not able to fill all of the missing data.\n',...
        'There are still %i days of missing data.\n'],sum(isnan(precip)));
end


end
