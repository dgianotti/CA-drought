function d1km=lldistkm(latlon1,latlon2)
% format: d1km = lldistkm(latlon1,latlon2)
% Distance:
% d1km: distance in km based on Haversine formula
% (Haversine: http://en.wikipedia.org/wiki/Haversine_formula)
%
% --Inputs:
%   latlon1: latlon of origin point [lat lon] (must be 1 x 2 vector)
%   latlon2: latlon of destination point [lat lon] (Can be an N x 2 arrray)
%
% --Outputs:
%   d1km: distance calculated by Haversine formula
%

if (size(latlon1) ~= [1,2])
    error('The first input to lldistkm should be a [1 x 2] vector. Aborting!');
elseif (size(latlon2,2) ~= 2)
    error('The second input to lldistkm shoud be an [N x 2] array. Aborting!');
end
     
radius=6371;
lat1=latlon1(1)*pi/180;
lat2=latlon2(:,1)*pi/180;
lon1=latlon1(2)*pi/180;
lon2=latlon2(:,2)*pi/180;

deltaLat=lat2-lat1;
deltaLon=lon2-lon1;

a=sin((deltaLat)/2).^2 + cos(lat1).*cos(lat2) .* sin(deltaLon/2).^2;
c=2*atan2(sqrt(a),sqrt(1-a));
d1km=radius*c;    %Haversine distance


end