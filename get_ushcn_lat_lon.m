function [lat, lon] = get_ushcn_lat_lon(id)

% Make sure you have the ushcn-stations file
if ~exist('ushcn-stations.txt','file')
    urlwrite('http://cdiac.ornl.gov/ftp/ushcn_daily/ushcn-stations.txt','ushcn-stations.txt');
end

% Load the ushcn file:
fmt = '%6s %1*s %8s %1*s %9s %65*s';
fid = fopen('ushcn-stations.txt');
F = textscan(fid, fmt, 'Delimiter','', 'HeaderLines',0,'Whitespace','');
fclose(fid);

LAT = str2double(F{2});
LON = str2double(F{3});

lat = LAT(strcmp(id,F{1}));
lon = LON(strcmp(id,F{1}));
end