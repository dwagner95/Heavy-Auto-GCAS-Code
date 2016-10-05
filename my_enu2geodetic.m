function [lat, lon, h] = my_enu2geodetic( ...
    xEast, yNorth, zUp, lat0, lon0, h0)
%GEODETIC2ENU Summary of this function goes here
%   Detailed explanation goes here

a = 6378137;
f = 1/298.257223563;

sinphi = sind(lat0);
cosphi = cosd(lat0);

e2 = f * (2 - f);
N  = a ./ sqrt(1 - e2 * sinphi.^2);
rho = (N + h0) .* cosphi;
ecefZ0 = (N*(1 - e2) + h0) .* sinphi;

ecefX0 = rho .* cosd(lon0);
ecefY0 = rho .* sind(lon0);

cosPhi = cosd(lat0);
sinPhi = sind(lat0);
cosLambda = cosd(lon0);
sinLambda = sind(lon0);

t = cosPhi .* zUp - sinPhi .* yNorth;
dz = sinPhi .* zUp + cosPhi .* yNorth;

dx = cosLambda .* t - sinLambda .* xEast;
dy = sinLambda .* t + cosLambda .* xEast;

x = ecefX0 + dx;
y = ecefY0 + dy;
z = ecefZ0 + dz;

rho = hypot(x,y);


b = (1 - f) * a;       % Semiminor axis
e2 = f * (2 - f);      % Square of (first) eccentricity
ep2 = e2 / (1 - e2);   % Square of second eccentricity

% Bowring's formula for initial parametric (beta) and geodetic
% (phi) latitudes
beta = atan2d(z, (1 - f) * rho);
lat = atan2d(z   + b * ep2 * sind(beta).^3,...
    rho - a * e2  * cosd(beta).^3);

% Fixed-point iteration with Bowring's formula
% (typically converges within two or three iterations)
betaNew = atan2d((1 - f)*sind(lat), cosd(lat));
count = 0;
while any(beta(:) ~= betaNew(:)) && count < 5
    beta = betaNew;
    lat = atan2d(z   + b * ep2 * sind(beta).^3,...
        rho - a * e2  * cosd(beta).^3);
    betaNew = atan2d((1 - f)*sind(lat), cosd(lat));
    count = count + 1;
end

% Ellipsoidal height from final value for latitude
sinphi = sind(lat);
N = a ./ sqrt(1 - e2 * sinphi.^2);
h = rho .* cosd(lat) + (z + e2 * N .* sinphi) .* sinphi - N;

lon = atan2d(y,x);


end

