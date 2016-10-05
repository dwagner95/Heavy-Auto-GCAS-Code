function [xEast, yNorth, zUp] = my_geodetic2enu( ...
    lat, lon, h, lat0, lon0, h0)
%GEODETIC2ENU Summary of this function goes here
%   Detailed explanation goes here

Eccentricity = 0.081819190842621;
a = 6378137;

e2 = Eccentricity^2;

s1 = sind(lat0);
c1 = cosd(lat0);

s2 = sind(lat);
c2 = cosd(lat);

p1 = c1 .* cosd(lon0);
p2 = c2 .* cosd(lon);

q1 = c1 .* sind(lon0);
q2 = c2 .* sind(lon);

w1 = 1 ./ sqrt(1 - e2 * s1.^2);
w2 = 1 ./ sqrt(1 - e2 * s2.^2);

deltaX =            a * (p2 .* w2 - p1 .* w1) + (h .* p2 - h0 .* p1);
deltaY =            a * (q2 .* w2 - q1 .* w1) + (h .* q2 - h0 .* q1);
deltaZ = (1 - e2) * a * (s2 .* w2 - s1 .* w1) + (h .* s2 - h0 .* s1);

cosPhi = cosd(lat0);
sinPhi = sind(lat0);
cosLambda = cosd(lon0);
sinLambda = sind(lon0);

t     =  cosLambda .* deltaX + sinLambda .* deltaY;
xEast = -sinLambda .* deltaX + cosLambda .* deltaY;

zUp    =  cosPhi .* t + sinPhi .* deltaZ;
yNorth = -sinPhi .* t + cosPhi .* deltaZ;

end

