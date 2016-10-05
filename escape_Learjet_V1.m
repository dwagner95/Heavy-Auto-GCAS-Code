function [Nz_vector_cmd, bank_vector_cmd, take_over ] = escape_learjet_v1(lat_pos,lon_pos,tf,jet_head,jet_gamma,alt,bubble )

% Maj Trombetta May 2015
clc;
%The primary use of this code will be the first version of the Learjet
%function that will be flown for 15A HAVE ESCAPE TMP.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%User Inputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%lat_star/lon_start/lat_end/lon_end: Degrees, I will more than likely not
%need these since we will just use current aircraft position and not the
%path created by these lat/lon



%This code will run the 5 avoidance path algorithm with ODE3.  It also
%takes into account a 0.5 second look-ahead buffer which is also propagated
%in the control vectors.  This code will work for a slow and medium speed
%aircraft (what we plan to flight test).  For a fast speed aircraft the 
%control and ODE functions will need to be changed.
%%%%%%%%%%%%%%%%%%%%%%Choose Flight Terrain Box%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These latitude and longitude should be set to encompass the entire
% flyable area, more than likely all of R-2508 or just R-2515.
latlim_S=35.1;
latlim_N=35.3;
lonlim_E=-117.3;
lonlim_W=-117.5;
%%%%%%%%%%%%%%%%%%%%%%Aircraft Flight Path%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E=referenceEllipsoid('earth'); %Sets Earth is the reference Ellipsoid
%%%%%%%%%%%%%%%%%Calculate Heading of Aircraft Flight%%%%%%%%%%%%%%%%%%%%%%
%head_calc=azimuth('rh',lat_start,lon_start,lat_end,lon_end,E);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%Set Initial Bubble (Avoidance) Condition and Time Window%%%%%%%%%%%
t0 = 0;                         %Iteration Start Time
%tf = 31;                        %Set to desired Path Propagation Time
look_ahead = 0.5;               %Look ahead buffer time (If You change this
% from 0.5 seconds, the control must also be changed to reflect the correct
% time.
%bubble =  300;                  %Set this in feet (Bubbles < 300 ft will
%use DTED2
%%%%%%%%%%%%%%%%%%%Enter Current Aircraft Parameters%%%%%%%%%%%%%%%%%%%%%%%
V=523.22;                       %Constant Velocity, ft/s (310 kts)
lat_pos = deg2rad(lat_pos);    %Initial Condition in Latitude (degrees)
lon_pos = deg2rad(lon_pos);    %Initial Condition in Longitude (degrees)
alt = ft2m(alt);               %MSL Altitude as read from Altimeter (ft)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%No User Inputs Below This Line%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%Build Terrain Map From DTED%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Map Read: GCAS Mtn in R-2508 at Edwards AFB CA
% DTED: SRT1f
latlim = [latlim_S latlim_N];
lonlim = [lonlim_W lonlim_E];
% Chooses higher fidelity DTED if bubble is small to ensure DTED Post
% Capture
if bubble < 300
    [Z, refvec, lat, lon] = dted('n35.dt2', 1, latlim, lonlim); % Elevation (m)
else
    [Z, refvec, lat, lon] = dted('n35.dt1', 1, latlim, lonlim); % Elevation (m)
end
%Z(Z == 0) = -1;
cellsperdegree = refvec(1);
step = 1/cellsperdegree;
Zlength = length(Z);            %meters
latN = refvec(2);     %pulls the northern lat position from refvec
lonW = refvec(3);     %pulls the Western lon position from refvec      
lonE = (Zlength-1)*step+lonW;  %calculates the opposite corner
latS = (Zlength-1)*(-step)+latN; %calculates the opposite corner

Xvector = lonW:step:lonE;
XLongMatrix = repmat(Xvector,Zlength,1); %deg

Yvector = latS:step:latN;
YLatMatrix = repmat(Yvector,Zlength,1)'; %deg
% load('gcasmap');
% load('blue_sp');
% load('gray_sp');
% load('green_sp');
% load('orange_sp');
% load('purple_sp');

%%% All Controls
load('bank_angle_lft');
load('bank_angle_rgt');
load('bank_angle_lft_up');
load('bank_angle_rgt_up');
load('Nz_pitch_up');
load('Nz_ba');
load('Nz_pitch');
%%% End Controls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%Polyfit for Length of a Degree of Longitude%%%%%%%%%%%%%%%%%%%
%Deg latitude
x_poly=0:10:90;
%Length of a degree of Geodetic Longitude at a specific latitude
y_poly=[111.32,109.64,104.65,96.49,85.39,71.70,55.80,38.19,19.39,0];
%order of polynomial fit
n_poly=5;
[P]=polyfit(x_poly,y_poly,n_poly);
%http://www.ncgia.ucsb.edu/education/curricula/giscc/units/u014/tables/table02.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Geoid Defined%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Average Geoid Height is -31.783 m averaged over entire surface plot. 
%HAE=MSL+G_Height, 'geoidHeight' MATLAB command
G_Height = -31.783;             % meters
alt = alt + G_Height;           %Current WGS-84 HAE A/C Altitude (m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%Set-up GriddedInterpolant%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xgrid = XLongMatrix'; %Sets up ndgrid format for the interpolate function
Ygrid = YLatMatrix';
Zgrid = (Z + G_Height)'; %Turns all DTED altitude values into HAE
F = griddedInterpolant(Xgrid,Ygrid,Zgrid,'spline'); %HAE Altitudes

%% Path Solutions: ODE45 for point mass 3-DOF Model
%Program to run aircraft equations of motion for 3-DOF point mass
%model and determine if a collision will occur and what path must be
%executed to avoid a collision.

%%%%%%%%%%Calculates required values for Sphere Analysis%%%%%%%%%%%%%%%%%%
bubble_meters = ft2m(bubble);%Puts bubble in meters for consistency
R2 = bubble_meters^2;        %Squared once for speed (m^2)
M = diag(1./[R2 R2 R2]);     %For transformation matrix (m^2)
km2deg = (1/111)*(1/3280.84);   %Conversion to change ft2deg in latitude
exitval = 0;
kkount = 0;
t_grnd_col_spy = zeros(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%calculates the HAE ground height @ current A/C Lat and Lon (griddedInterpolant)
h0=F(rad2deg(lon_pos),rad2deg(lat_pos)); % F(Longitude,Latitude)

% Find current position in ENU Reference Frame
[x_init,y_init,z_init]=geodetic2enu(...
    lat_pos,lon_pos,alt,...
    lat_pos,lon_pos,h0,...
    E,'radians'); %Answer in meters

% Initialize A/C state values for Propagation with ODE45
x_init = m2ft(x_init);     % Initial x position, ft
y_init = m2ft(y_init);     % Initial y position, ft
z_init = m2ft(z_init);     % Initial z position, ft (Altitude)
gam_init =  jet_gamma;     % Initial flight path angle, degrees
chi_init = compass2pol(jet_head);    % Initial heading angle

%Initial conditions column vector [x;y;z;gam;chi] [ft,ft,ft,deg,deg]
Init_con= [x_init;y_init;z_init;deg2rad(gam_init);deg2rad(chi_init)]; 

%Runs ODE45 for to find the A/C state for each maneuver
tspan=t0:.3:(tf+look_ahead);
% Adding .5 alots for the look-ahead buffer
[state_fwd]=ode3(@ac_state_fwd,tspan,Init_con,Nz_pitch,V);
[state_lft]=ode3(@ac_state_lft,tspan,Init_con,Nz_ba,bank_angle_lft,V);
[state_rgt]=ode3(@ac_state_rgt,tspan,Init_con,Nz_ba,bank_angle_rgt,V);
[state_lft_up]=ode3(@ac_state_lft_up,tspan,Init_con,Nz_pitch_up,bank_angle_lft_up,V);
[state_rgt_up]=ode3(@ac_state_rgt_up,tspan,Init_con,Nz_pitch_up,bank_angle_rgt_up,V);

%Define t_fwd, t_rgt, t_lft so that I don't have to change all the logic,
%they are unnecessary with a predefined timespan but it's easier to
%redifine them to that timespan in lieu of changing all the code.
t_fwd = tspan;
t_lft = tspan;
t_rgt = tspan;
t_lft_up = tspan;
t_rgt_up = tspan;

%Converts state solutions back to degrees
state_fwd=[state_fwd(:,1),state_fwd(:,2),state_fwd(:,3),rad2deg(state_fwd(:,4)),...
    rad2deg(state_fwd(:,5))];

state_lft=[state_lft(:,1),state_lft(:,2),state_lft(:,3),rad2deg(state_lft(:,4)),...
    rad2deg(state_lft(:,5))];

state_rgt=[state_rgt(:,1),state_rgt(:,2),state_rgt(:,3),rad2deg(state_rgt(:,4)),...
    rad2deg(state_rgt(:,5))];

state_lft_up=[state_lft_up(:,1),state_lft_up(:,2),state_lft_up(:,3),rad2deg(state_lft_up(:,4)),...
    rad2deg(state_lft_up(:,5))];

state_rgt_up=[state_rgt_up(:,1),state_rgt_up(:,2),state_rgt_up(:,3),rad2deg(state_rgt_up(:,4)),...
    rad2deg(state_rgt_up(:,5))];

%Convert propagated position back into Lat and Lon (Fwd Solution)
[lat_new_fwd,lon_new_fwd,alt_new_fwd]=enu2geodetic(ft2m(state_fwd(:,1)),...
    ft2m(state_fwd(:,2)),ft2m(state_fwd(:,3)),lat_pos, lon_pos,h0,E,...
    'radians');
lat_new_fwd=rad2deg(lat_new_fwd);
lon_new_fwd=rad2deg(lon_new_fwd);

%Convert propagated position back into Lat and Lon (Lft Solution)
[lat_new_lft,lon_new_lft,alt_new_lft]=enu2geodetic(ft2m(state_lft(:,1)),...
    ft2m(state_lft(:,2)),ft2m(state_lft(:,3)),lat_pos, lon_pos,h0,E,...
    'radians');
lat_new_lft=rad2deg(lat_new_lft);
lon_new_lft=rad2deg(lon_new_lft);

%Convert propagated position back into Lat and Lon (Rgt Solution)
[lat_new_rgt,lon_new_rgt,alt_new_rgt]=enu2geodetic(ft2m(state_rgt(:,1)),...
    ft2m(state_rgt(:,2)),ft2m(state_rgt(:,3)),lat_pos, lon_pos,h0,E,...
    'radians');
lat_new_rgt=rad2deg(lat_new_rgt);
lon_new_rgt=rad2deg(lon_new_rgt);

%Convert propagated position back into Lat and Lon (Lft Up Solution)
[lat_new_lft_up,lon_new_lft_up,alt_new_lft_up]=enu2geodetic(ft2m(state_lft_up(:,1)),...
    ft2m(state_lft_up(:,2)),ft2m(state_lft_up(:,3)),lat_pos, lon_pos,h0,E,...
    'radians');
lat_new_lft_up=rad2deg(lat_new_lft_up);
lon_new_lft_up=rad2deg(lon_new_lft_up);

%Convert propagated position back into Lat and Lon (Rgt Up Solution)
[lat_new_rgt_up,lon_new_rgt_up,alt_new_rgt_up]=enu2geodetic(ft2m(state_rgt_up(:,1)),...
    ft2m(state_rgt_up(:,2)),ft2m(state_rgt_up(:,3)),lat_pos, lon_pos,h0,E,...
    'radians');
lat_new_rgt_up=rad2deg(lat_new_rgt_up);
lon_new_rgt_up=rad2deg(lon_new_rgt_up);

% Logic to check the fwd solution for terrain collision
% (Fwd Solution)
sp_dist_fwd = zeros(1,length(t_fwd));
sp_dist = zeros(1,1);
lat_lon_alt = zeros(1,1);
for i=1:length(t_fwd)
    
%Finds the distance between a deg of longitude at a specific Lat in km
long_width = polyval(P,lat_new_fwd(i));
%Conversion from km to degrees for longitude distance
lon_bubble = inv(long_width)*inv(3280.84)*bubble; %deg Long at Lat position
    
%finds the column indices Longitude DTED Posts
long_find=find(Xvector <= (lon_new_fwd(i)+(lon_bubble))...
& Xvector >= (lon_new_fwd(i)-(lon_bubble)));

%Finds the row indices Latitude DTED Posts
lat_find=find(Yvector <= (lat_new_fwd(i)+(km2deg*bubble))...
& Yvector >= (lat_new_fwd(i)-(km2deg*bubble)));

%Breaks current loop iteration if no DTED Posts found
if isempty(long_find)==1 || isempty(lat_find)==1
    sp_dist_fwd(i) = 8888;
    continue
end

%Finds Long/Lat values from complete Long/Lat Matrices
long_mat_find = XLongMatrix(lat_find,long_find);
lat_mat_find = YLatMatrix(lat_find,long_find);
Z_mat_find = Zgrid(long_find,lat_find); %(long,lat) because Zgrid is transposed

%Build lat_lon_alt DTED posts for this timestep
kount = numel(long_mat_find);
for k=1:kount
lat_lon_alt(k,1:3)=[lat_mat_find(k),long_mat_find(k),Z_mat_find(k)];
%vpa(lat_lon_alt)
end

%convert DTED to ENU
[x_fwd,y_fwd,z_fwd]=geodetic2enu(deg2rad(lat_lon_alt(:,1)),deg2rad(lat_lon_alt(:,2)),...
    lat_lon_alt(:,3),lat_pos, lon_pos,h0,E,'radians'); %Answer in meters

%Build matrices required to check if DTED are inside bubble
xyz_dted = [x_fwd,y_fwd,z_fwd];
AC_pos = [ft2m(state_fwd(i,1)),ft2m(state_fwd(i,2)),ft2m(state_fwd(i,3))];

%Checks if DTED Post height is inside sphere XMX'
for j=1:length(x_fwd)
   sp_dist(j) =  (xyz_dted(j,:)- AC_pos)*M*(xyz_dted(j,:)- AC_pos)'; 
end

%Finds the smallest distance for a given timestep
sp_dist_fwd(i) = min(sp_dist);
clear sp_dist;
clear lat_lon_alt;
end

%Calculates the fwd time of a collision occurring.
t_grnd_col_fwd = t_fwd(sp_dist_fwd <=1);

if isempty(t_grnd_col_fwd) == 0 %Meaning: not empty
    t_grnd_col_fwd = t_grnd_col_fwd(1);
else
    %Fills an empty matrix with -1 to delineate that no collision occured
    t_grnd_col_fwd = -1;
end

% Logic to check the lft solution for terrain collision
% (Lft Solution)
sp_dist_lft = zeros(1,length(t_lft));
sp_dist = zeros(1,1);
lat_lon_alt = zeros(1,1);
for i=1:length(t_lft)
    
%Finds the distance between a deg of longitude at a specific Lat in km
long_width = polyval(P,lat_new_lft(i));
%Conversion from km to degrees for longitude distance
lon_bubble = inv(long_width)*inv(3280.84)*bubble; %deg Long at Lat position
    
    %finds the column indices Longitude DTED Posts
long_find=find(Xvector <= (lon_new_lft(i)+(lon_bubble))...
& Xvector >= (lon_new_lft(i)-(lon_bubble)));

%Finds the row indices Latitude DTED Posts
lat_find=find(Yvector <= (lat_new_lft(i)+(km2deg*bubble))...
& Yvector >= (lat_new_lft(i)-(km2deg*bubble)));

%Breaks current loop iteration if no DTED Posts found
if isempty(long_find)==1 || isempty(lat_find)==1
    sp_dist_lft(i) = 8888;
    continue
end

%Finds Long/Lat values from complete Long/Lat Matrices
long_mat_find = XLongMatrix(lat_find,long_find);
lat_mat_find = YLatMatrix(lat_find,long_find);
Z_mat_find = Zgrid(long_find,lat_find); %(long,lat) because Zgrid is transposed

%Build lat_lon_alt DTED posts for this timestep
kount = numel(long_mat_find);
for k=1:kount
lat_lon_alt(k,1:3)=[lat_mat_find(k),long_mat_find(k),Z_mat_find(k)];
end

%convert DTED to ENU
[x_lft,y_lft,z_lft]=geodetic2enu(deg2rad(lat_lon_alt(:,1)),deg2rad(lat_lon_alt(:,2)),...
    lat_lon_alt(:,3),lat_pos, lon_pos,h0,E,'radians'); %Answer in meters

%Build matrices required to check if DTED are inside bubble
xyz_dted = [x_lft,y_lft,z_lft];
AC_pos = [ft2m(state_lft(i,1)),ft2m(state_lft(i,2)),ft2m(state_lft(i,3))];

%Checks if DTED Post height is inside sphere XMX'
for j=1:length(x_lft)
   sp_dist(j) =  (xyz_dted(j,:)- AC_pos)*M*(xyz_dted(j,:)- AC_pos)';
end

%Finds the smallest distance for a given timestep
sp_dist_lft(i) = min(sp_dist);
clear sp_dist;
clear lat_lon_alt;
end

%Calculates the lft time of a collision occurring.
t_grnd_col_lft = t_lft(sp_dist_lft <= 1);

if isempty(t_grnd_col_lft) == 0 %Meaning: not empty
    t_grnd_col_lft = t_grnd_col_lft(1);
else
    %Fills an empty matrix with -1 to delineate that no collision occured
    t_grnd_col_lft = -1;
end

% Logic to check the rgt solution for terrain collision
% (Rgt Solution)
sp_dist_rgt = zeros(1,length(t_rgt));
sp_dist = zeros(1,1);
lat_lon_alt = zeros(1,1);
for i=1:length(t_rgt)
    
%Finds the distance between a deg of longitude at a specific Lat in km
long_width = polyval(P,lat_new_rgt(i));
%Conversion from km to degrees for longitude distance
lon_bubble = inv(long_width)*inv(3280.84)*bubble; %deg Long at Lat position
    
%finds the column indices Longitude DTED Posts
long_find=find(Xvector <= (lon_new_rgt(i)+(lon_bubble))...
& Xvector >= (lon_new_rgt(i)-(lon_bubble)));

%Finds the row indices Latitude DTED Posts
lat_find=find(Yvector <= (lat_new_rgt(i)+(km2deg*bubble))...
& Yvector >= (lat_new_rgt(i)-(km2deg*bubble)));

%Breaks current loop iteration if no DTED Posts found
if isempty(long_find)==1 || isempty(lat_find)==1
    sp_dist_rgt(i) = 8888;
    continue
end

%Finds Long/Lat values from complete Long/Lat Matrices
long_mat_find = XLongMatrix(lat_find,long_find);
lat_mat_find = YLatMatrix(lat_find,long_find);
Z_mat_find = Zgrid(long_find,lat_find); %(long,lat) because Zgrid is transposed

%Build lat_lon_alt DTED posts for this timestep
kount = numel(long_mat_find);
for k=1:kount
lat_lon_alt(k,1:3)=[lat_mat_find(k),long_mat_find(k),Z_mat_find(k)];
end

%convert DTED to ENU
[x_rgt,y_rgt,z_rgt]=geodetic2enu(deg2rad(lat_lon_alt(:,1)),deg2rad(lat_lon_alt(:,2)),...
    lat_lon_alt(:,3),lat_pos, lon_pos,h0,E,'radians'); %Answer in meters

%Build matrices required to check if DTED are inside bubble
xyz_dted = [x_rgt,y_rgt,z_rgt];
AC_pos = [ft2m(state_rgt(i,1)),ft2m(state_rgt(i,2)),ft2m(state_rgt(i,3))];

%Checks if DTED Post height is inside sphere XMX'
for j=1:length(x_rgt)
   sp_dist(j) =  (xyz_dted(j,:)- AC_pos)*M*(xyz_dted(j,:)- AC_pos)'; 
end

%Finds the smallest distance for a given timestep
sp_dist_rgt(i) = min(sp_dist);
clear sp_dist;
clear lat_lon_alt;
end

%Calculates the rgt time of a collision occurring.
t_grnd_col_rgt = t_rgt(sp_dist_rgt <= 1);

if isempty(t_grnd_col_rgt) == 0 %Meaning: not empty
    t_grnd_col_rgt = t_grnd_col_rgt(1);
else
    %Fills an empty matrix with -1 to delineate that no collision occured
    t_grnd_col_rgt = -1;
end

%%
% Logic to check the lft_up solution for terrain collision
% (lft_up Solution)
sp_dist_lft_up = zeros(1,length(t_lft_up));
sp_dist = zeros(1,1);
lat_lon_alt = zeros(1,1);
for i=1:length(t_lft_up)
    
%Finds the distance between a deg of longitude at a specific Lat in km
long_width = polyval(P,lat_new_lft_up(i));
%Conversion from km to degrees for longitude distance
lon_bubble = inv(long_width)*inv(3280.84)*bubble; %deg Long at Lat position
    
%finds the column indices Longitude DTED Posts
long_find=find(Xvector <= (lon_new_lft_up(i)+(lon_bubble))...
& Xvector >= (lon_new_lft_up(i)-(lon_bubble)));

%Finds the row indices Latitude DTED Posts
lat_find=find(Yvector <= (lat_new_lft_up(i)+(km2deg*bubble))...
& Yvector >= (lat_new_lft_up(i)-(km2deg*bubble)));

%Breaks current loop iteration if no DTED Posts found
if isempty(long_find)==1 || isempty(lat_find)==1
    sp_dist_fwd(i) = 8888;
    continue
end

%Finds Long/Lat values from complete Long/Lat Matrices
long_mat_find = XLongMatrix(lat_find,long_find);
lat_mat_find = YLatMatrix(lat_find,long_find);
Z_mat_find = Zgrid(long_find,lat_find); %(long,lat) because Zgrid is transposed

%Build lat_lon_alt DTED posts for this timestep
kount = numel(long_mat_find);
for k=1:kount
lat_lon_alt(k,1:3)=[lat_mat_find(k),long_mat_find(k),Z_mat_find(k)];
%vpa(lat_lon_alt)
end

%convert DTED to ENU
[x_lft_up,y_lft_up,z_lft_up]=geodetic2enu(deg2rad(lat_lon_alt(:,1)),deg2rad(lat_lon_alt(:,2)),...
    lat_lon_alt(:,3),lat_pos, lon_pos,h0,E,'radians'); %Answer in meters

%Build matrices required to check if DTED are inside bubble
xyz_dted = [x_lft_up,y_lft_up,z_lft_up];
AC_pos = [ft2m(state_lft_up(i,1)),ft2m(state_lft_up(i,2)),ft2m(state_lft_up(i,3))];

%Checks if DTED Post height is inside sphere XMX'
for j=1:length(x_lft_up)
   sp_dist(j) =  (xyz_dted(j,:)- AC_pos)*M*(xyz_dted(j,:)- AC_pos)'; 
end

%Finds the smallest distance for a given timestep
sp_dist_lft_up(i) = min(sp_dist);
clear sp_dist;
clear lat_lon_alt;
end

%Calculates the fwd time of a collision occurring.
t_grnd_col_lft_up = t_fwd(sp_dist_lft_up <=1);

if isempty(t_grnd_col_lft_up) == 0 %Meaning: not empty
    t_grnd_col_lft_up = t_grnd_col_lft_up(1);
else
    %Fills an empty matrix with -1 to delineate that no collision occured
    t_grnd_col_lft_up = -1;
end

%%
% Logic to check the rgt_up solution for terrain collision
% (rgt_up Solution)
sp_dist_rgt_up = zeros(1,length(t_rgt_up));
sp_dist = zeros(1,1);
lat_lon_alt = zeros(1,1);
for i=1:length(t_rgt_up)
    
%Finds the distance between a deg of longitude at a specific Lat in km
long_width = polyval(P,lat_new_rgt_up(i));
%Conversion from km to degrees for longitude distance
lon_bubble = inv(long_width)*inv(3280.84)*bubble; %deg Long at Lat position
    
%finds the column indices Longitude DTED Posts
long_find=find(Xvector <= (lon_new_rgt_up(i)+(lon_bubble))...
& Xvector >= (lon_new_rgt_up(i)-(lon_bubble)));

%Finds the row indices Latitude DTED Posts
lat_find=find(Yvector <= (lat_new_rgt_up(i)+(km2deg*bubble))...
& Yvector >= (lat_new_rgt_up(i)-(km2deg*bubble)));

%Breaks current loop iteration if no DTED Posts found
if isempty(long_find)==1 || isempty(lat_find)==1
    sp_dist_fwd(i) = 8888;
    continue
end

%Finds Long/Lat values from complete Long/Lat Matrices
long_mat_find = XLongMatrix(lat_find,long_find);
lat_mat_find = YLatMatrix(lat_find,long_find);
Z_mat_find = Zgrid(long_find,lat_find); %(long,lat) because Zgrid is transposed

%Build lat_lon_alt DTED posts for this timestep
kount = numel(long_mat_find);
for k=1:kount
lat_lon_alt(k,1:3)=[lat_mat_find(k),long_mat_find(k),Z_mat_find(k)];
%vpa(lat_lon_alt)
end

%convert DTED to ENU
[x_rgt_up,y_rgt_up,z_rgt_up]=geodetic2enu(deg2rad(lat_lon_alt(:,1)),deg2rad(lat_lon_alt(:,2)),...
    lat_lon_alt(:,3),lat_pos, lon_pos,h0,E,'radians'); %Answer in meters

%Build matrices required to check if DTED are inside bubble
xyz_dted = [x_rgt_up,y_rgt_up,z_rgt_up];
AC_pos = [ft2m(state_rgt_up(i,1)),ft2m(state_rgt_up(i,2)),ft2m(state_rgt_up(i,3))];

%Checks if DTED Post height is inside sphere XMX'
for j=1:length(x_rgt_up)
   sp_dist(j) =  (xyz_dted(j,:)- AC_pos)*M*(xyz_dted(j,:)- AC_pos)'; 
end

%Finds the smallest distance for a given timestep
sp_dist_rgt_up(i) = min(sp_dist);
clear sp_dist;
clear lat_lon_alt;
end

%Calculates the fwd time of a collision occurring.
t_grnd_col_rgt_up = t_fwd(sp_dist_rgt_up <=1);

if isempty(t_grnd_col_rgt_up) == 0 %Meaning: not empty
    t_grnd_col_rgt_up = t_grnd_col_rgt_up(1);
else
    %Fills an empty matrix with -1 to delineate that no collision occured
    t_grnd_col_rgt_up = -1;
end


%Builds time vector of potential collisions [fwd,lft,rgt]
t_grnd_col_all = [t_grnd_col_fwd,t_grnd_col_lft,t_grnd_col_rgt,t_grnd_col_lft_up,t_grnd_col_rgt_up];

%%%%%%%%%%%%%%%%%%%%SPY Plot (Figure 2) Command Code%%%%%%%%%%%%%%%%%%%%%%%
t_grnd_col_spy(c,1:5)=t_grnd_col_all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Converts Altitudes back to MSL from HAE (meters) MSL = HAE - G_Height
alt_new_fwd = alt_new_fwd - G_Height;
alt_new_lft = alt_new_lft - G_Height;
alt_new_rgt = alt_new_rgt - G_Height;
alt_new_lft_up = alt_new_lft_up - G_Height;
alt_new_rgt_up = alt_new_rgt_up - G_Height;


%Need to find a way to send the command vectors to the learjet.
% Calculates appropriate path and when to execute it.
if all(t_grnd_col_all(:) >= 0)  
    t_execute = max(t_grnd_col_all);
    if t_grnd_col_fwd >= t_grnd_col_lft && t_grnd_col_fwd >= t_grnd_col_rgt...
            && t_grnd_col_fwd >= t_grnd_col_lft_up && t_grnd_col_fwd >= t_grnd_col_rgt_up
           fprintf('Execute Forward Path \n')
           Nz_vector_cmd = Nz_pitch; %Execute Forward path (Nz_pitch.mat)
           flagg=1; %Breaks all the loops and ends the program
           break
       elseif t_grnd_col_lft > t_grnd_col_fwd && t_grnd_col_lft > t_grnd_col_rgt...
               && t_grnd_col_lft > t_grnd_col_lft_up && t_grnd_col_lft > t_grnd_col_rgt_up
           fprintf('Execute Left Path \n')
           Nz_vector_cmd = Nz_ba; %Execute Left Path (Nz_ba.mat, bank_angle_lft.mat)
           bank_vector_cmd = bank_angle_lft;
           flagg=1; %Breaks all the loops and ends the program
           break
       elseif t_grnd_col_rgt > t_grnd_col_fwd && t_grnd_col_rgt > t_grnd_col_lft...
               && t_grnd_col_rgt > t_grnd_col_lft_up && t_grnd_col_rgt > t_grnd_col_rgt_up
           fprintf('Execute Right Path \n')
           Nz_vector_cmd = Nz_ba; %Execute Right Path (Nz_ba.mat, bank_angle_rgt.mat)
           bank_vector_cmd = bank_angle_rgt;
           flagg=1; %Breaks all the loops and ends the program
           break
       elseif t_grnd_col_lft_up > t_grnd_col_fwd && t_grnd_col_lft_up >= t_grnd_col_rgt...
               && t_grnd_col_lft_up >= t_grnd_col_lft && t_grnd_col_lft_up >= t_grnd_col_rgt_up
           fprintf('Execute Left-Up Path \n')
           Nz_vector_cmd = Nz_pitch_up;%Execute Left-Up Path (Nz_pitch_up.mat, bank_angle_lft_up.mat)
           bank_vector_cmd = bank_angle_lft_up;
           flagg=1; %Breaks all the loops and ends the program 
           break
       elseif t_grnd_col_rgt_up > t_grnd_col_fwd && t_grnd_col_rgt_up >= t_grnd_col_rgt...
               && t_grnd_col_rgt_up >= t_grnd_col_lft && t_grnd_col_rgt_up > t_grnd_col_lft_up
           fprintf('Execute Right-Up Path \n')
           Nz_vector_cmd = Nz_pitch_up;%Execute Right-Up Path (Nz_pitch_up.mat, bank_angle_rgt_up.mat)
           bank_vector_cmd = bank_angle_rgt_up;
           flagg=1; %Breaks all the loops and ends the program
           break
    end   
else 
   fprintf('No Automated Path Deviation Required \n') 
end

if(flagg==1)
    take_over = 1; %flag that sets control vector to be sent to Learjet
    break
end
   
