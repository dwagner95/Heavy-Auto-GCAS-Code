% Maj Trombetta November 2014 
close all; clc;
load('LowAlt_PilotPath_3200ftLevel_Shorter2_LLA.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%User Inputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code will run the 5 avoidance path algorithm with ODE3.  It also
%takes into account a 0.5 second look-ahead buffer which is also propagated
%in the control vectors.  This code will work for a slow and medium speed
%aircraft (what we plan to flight test).  For a fast speed aircraft the 
%control and ODE functions will need to be changed.
%%%%%%%%%%%%%%%%%%%%%%Choose Flight Terrain Box%%%%%%%%%%%%%%%%%%%%%%%%%%%%
latlim_S=35.1;
latlim_N=35.3;
lonlim_E=-117.3;
lonlim_W=-117.5;
%%%%%%%%%%%%%%%%%%%%%%Aircraft Flight Path%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sp_step  = 5;           %Steps size of sphere plots (bigger = faster build)
N = 30;%length(xVec_Lon);                  %Number of iteration steps
E=referenceEllipsoid('earth'); %Sets Earth is the reference Ellipsoid
%%%%%%%%%%%%%%%%%Calculate Heading of Aircraft Flight%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%Set Initial Bubble (Avoidance) Condition and Time Window%%%%%%%%%%%
t0 = 0;                       %Iteration Start Time
tf = 45;                      %Set to desired Path Propagation Time
look_ahead = 0;               %Look ahead buffer time (If You change this
% from 0.5 seconds, the control must also be changed to reflect the correct
% time.
bubble =  350;                  %Set this in feet (Bubbles < 300 ft will
%use DTED2
%%%%%%%%%%%%%%%%%%%Enter Current Aircraft Parameters%%%%%%%%%%%%%%%%%%%%%%%
V=354.44;                       % Constant Velocity, ft/s (210 kts)
lat_pos = deg2rad(yVec_Lat);    %Initial Condition in Latitude (degrees)
lon_pos = deg2rad(xVec_Lon);    %Initial Condition in Longitude (degrees)
jet_head = pol2compass(rad2deg(heading));%Input A/C heading (compass heading, deg)
jet_gamma = gamma;                  %Enter current Climb Angle (deg)
alt = Zac;               %MSL Altitude as read from Altimeter (ft)
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
%%% All Colormaps
load('gcasmap');
load('blue_sp');
load('gray_sp');
load('green_sp');
load('orange_sp');
load('purple_sp');
%%% End Colormaps

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
                                % MSL = HAE - G
                                % HAE = MSL + G

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%Set-up GriddedInterpolant%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xgrid = XLongMatrix'; %Sets up ndgrid format for the interpolate function
Ygrid = YLatMatrix';
Zgrid = (Z + G_Height)'; %Turns all DTED altitude values into HAE
F = griddedInterpolant(Xgrid,Ygrid,Zgrid,'spline'); %HAE Altitudes
% figure (3)
% surf(Xgrid,Ygrid,Zgrid)
%mesh(Xgrid,Ygrid,Zgrid)
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
while ~exitval
   
for c=1:N
kkount=kkount+1;
fprintf('Iteration %.0f \n',c);
tic
%calculates the HAE ground height @ current A/C Lat and Lon (griddedInterpolant)
h0=F(rad2deg(lon_pos(c)),rad2deg(lat_pos(c))); % F(Longitude,Latitude)

% Find current position in ENU Reference Frame
[x_init,y_init,z_init]=geodetic2enu(lat_pos(c),lon_pos(c),alt(c),lat_pos(c),lon_pos(c),...
    h0,E,'radians'); %Answer in meters

% Initialize A/C state values for Propagation with ODE45
x_init = m2ft(x_init);     % Initial x position, ft
y_init = m2ft(y_init);     % Initial y position, ft
z_init = m2ft(z_init);     % Initial z position, ft (Altitude)
gam_init =  jet_gamma(c);     % Initial flight path angle, degrees
chi_init = compass2pol(jet_head(c));    % Initial heading angle

%Initial conditions column vector [x;y;z;gam;chi] [ft,ft,ft,rad,rad]
Init_con= [x_init;y_init;z_init;gam_init;deg2rad(chi_init)]; 

%Runs ODE45 for to find the A/C state for each maneuver
tspan=t0:.5:(tf+look_ahead);
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
    ft2m(state_fwd(:,2)),ft2m(state_fwd(:,3)),lat_pos(c), lon_pos(c),h0,E,...
    'radians');
lat_new_fwd=rad2deg(lat_new_fwd);
lon_new_fwd=rad2deg(lon_new_fwd);

%Convert propagated position back into Lat and Lon (Lft Solution)
[lat_new_lft,lon_new_lft,alt_new_lft]=enu2geodetic(ft2m(state_lft(:,1)),...
    ft2m(state_lft(:,2)),ft2m(state_lft(:,3)),lat_pos(c), lon_pos(c),h0,E,...
    'radians');
lat_new_lft=rad2deg(lat_new_lft);
lon_new_lft=rad2deg(lon_new_lft);

%Convert propagated position back into Lat and Lon (Rgt Solution)
[lat_new_rgt,lon_new_rgt,alt_new_rgt]=enu2geodetic(ft2m(state_rgt(:,1)),...
    ft2m(state_rgt(:,2)),ft2m(state_rgt(:,3)),lat_pos(c), lon_pos(c),h0,E,...
    'radians');
lat_new_rgt=rad2deg(lat_new_rgt);
lon_new_rgt=rad2deg(lon_new_rgt);

%Convert propagated position back into Lat and Lon (Lft Up Solution)
[lat_new_lft_up,lon_new_lft_up,alt_new_lft_up]=enu2geodetic(ft2m(state_lft_up(:,1)),...
    ft2m(state_lft_up(:,2)),ft2m(state_lft_up(:,3)),lat_pos(c), lon_pos(c),h0,E,...
    'radians');
lat_new_lft_up=rad2deg(lat_new_lft_up);
lon_new_lft_up=rad2deg(lon_new_lft_up);

%Convert propagated position back into Lat and Lon (Rgt Up Solution)
[lat_new_rgt_up,lon_new_rgt_up,alt_new_rgt_up]=enu2geodetic(ft2m(state_rgt_up(:,1)),...
    ft2m(state_rgt_up(:,2)),ft2m(state_rgt_up(:,3)),lat_pos(c), lon_pos(c),h0,E,...
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
    lat_lon_alt(:,3),lat_pos(c), lon_pos(c),h0,E,'radians'); %Answer in meters

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
    lat_lon_alt(:,3),lat_pos(c), lon_pos(c),h0,E,'radians'); %Answer in meters

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
    lat_lon_alt(:,3),lat_pos(c), lon_pos(c),h0,E,'radians'); %Answer in meters

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
    lat_lon_alt(:,3),lat_pos(c), lon_pos(c),h0,E,'radians'); %Answer in meters

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
    lat_lon_alt(:,3),lat_pos(c), lon_pos(c),h0,E,'radians'); %Answer in meters

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

% Calculates appropriate path and when to execute it.
if all(t_grnd_col_all(:) >= 0)  
    t_execute = max(t_grnd_col_all);
    if t_grnd_col_fwd >= t_grnd_col_lft && t_grnd_col_fwd >= t_grnd_col_rgt...
            && t_grnd_col_fwd >= t_grnd_col_lft_up && t_grnd_col_fwd >= t_grnd_col_rgt_up
           fprintf('Execute Forward Path \n')
       elseif t_grnd_col_lft > t_grnd_col_fwd && t_grnd_col_lft > t_grnd_col_rgt...
               && t_grnd_col_lft > t_grnd_col_lft_up && t_grnd_col_lft > t_grnd_col_rgt_up
           fprintf('Execute Left Path \n')
       elseif t_grnd_col_rgt > t_grnd_col_fwd && t_grnd_col_rgt > t_grnd_col_lft...
               && t_grnd_col_rgt > t_grnd_col_lft_up && t_grnd_col_rgt > t_grnd_col_rgt_up
           fprintf('Execute Right Path \n')
       elseif t_grnd_col_lft_up > t_grnd_col_fwd && t_grnd_col_lft_up > t_grnd_col_rgt...
               && t_grnd_col_lft_up > t_grnd_col_lft && t_grnd_col_lft_up > t_grnd_col_rgt_up
           fprintf('Execute Left-Up Path \n')
       elseif t_grnd_col_rgt_up > t_grnd_col_fwd && t_grnd_col_rgt_up > t_grnd_col_rgt...
               && t_grnd_col_rgt_up > t_grnd_col_lft && t_grnd_col_rgt_up > t_grnd_col_lft_up
           fprintf('Execute Right-Up Path \n')
    end   
else 
   fprintf('No Automated Path Deviation Required \n') 
end

%Create Collision Report
fprintf('Collision Report: \n')
if all(t_grnd_col_all == -1)
    fprintf('No Path Collided with Terrain \n')
end  
if t_grnd_col_all(1) > -1
    fprintf('Forward Path Collided %.2f seconds from start \n',t_grnd_col_fwd);
end
if t_grnd_col_all(2) > -1
    fprintf('Left Path Collided %.2f seconds from start \n',t_grnd_col_lft);
end
if t_grnd_col_all(3) > -1
    fprintf('Right Path Collided %.2f seconds from start \n',t_grnd_col_rgt);
end
if t_grnd_col_all(4) > -1
    fprintf('Left-Up Collided %.2f seconds from start \n',t_grnd_col_lft_up);
end
if t_grnd_col_all(5) > -1
    fprintf('Right-Up Path Collided %.2f seconds from start \n',t_grnd_col_rgt_up);
end

%Converts Altitudes back to MSL from HAE (meters) MSL = HAE - G_Height
alt_new_fwd = alt_new_fwd - G_Height;
alt_new_lft = alt_new_lft - G_Height;
alt_new_rgt = alt_new_rgt - G_Height;
alt_new_lft_up = alt_new_lft_up - G_Height;
alt_new_rgt_up = alt_new_rgt_up - G_Height;
toc 
%%%%%%%%%%%%Surface Plots of the Terrain and Projected Path%%%%%%%%%%%%%%%%
%Creation of Spheres for plot
[x_sp,y_sp,z_sp]=sphere(10);
figure(1)
if kkount == 1
h1=surf(XLongMatrix,YLatMatrix,Z,'edgecolor','none');
colormap(gcasmap)
freezeColors
camlight(45,45)
hold on
end
%Forward Path Spheres
for i=1:sp_step:length(state_fwd)
    sphereE=ft2m((x_sp.*bubble) + state_fwd(i,1));
    sphereN=ft2m((y_sp.*bubble) + state_fwd(i,2));
    sphereU=ft2m((z_sp.*bubble) + state_fwd(i,3));
    
    [lat_sp, lon_sp, alt_sp] = enu2geodetic(sphereE,sphereN,sphereU,...
        lat_pos(c), lon_pos(c),h0,E,'radians');
    
    lat_sp=rad2deg(lat_sp);
    lon_sp=rad2deg(lon_sp);
    alt_sp = alt_sp - G_Height;
    h2=mesh(lon_sp,lat_sp,alt_sp);
    
end
    colormap(gray_sp)
    freezeColors
    
%Left Path Spheres
for i=1:sp_step:length(state_lft)
    sphereE=ft2m((x_sp.*bubble) + state_lft(i,1));
    sphereN=ft2m((y_sp.*bubble) + state_lft(i,2));
    sphereU=ft2m((z_sp.*bubble) + state_lft(i,3));
    
    [lat_sp, lon_sp, alt_sp] = enu2geodetic(sphereE,sphereN,sphereU,...
        lat_pos(c), lon_pos(c),h0,E,'radians');
    
    lat_sp=rad2deg(lat_sp);
    lon_sp=rad2deg(lon_sp);
    alt_sp = alt_sp - G_Height;
    h3=mesh(lon_sp,lat_sp,alt_sp);
    
end
    colormap(green_sp)
    freezeColors
    
%Left_Up Path Spheres
for i=1:sp_step:length(state_lft_up)
    sphereE=ft2m((x_sp.*bubble) + state_lft_up(i,1));
    sphereN=ft2m((y_sp.*bubble) + state_lft_up(i,2));
    sphereU=ft2m((z_sp.*bubble) + state_lft_up(i,3));
    
    [lat_sp, lon_sp, alt_sp] = enu2geodetic(sphereE,sphereN,sphereU,...
        lat_pos(c), lon_pos(c),h0,E,'radians');
    
    lat_sp=rad2deg(lat_sp);
    lon_sp=rad2deg(lon_sp);
    alt_sp = alt_sp - G_Height;
    h4=mesh(lon_sp,lat_sp,alt_sp);
    
end
    colormap(orange_sp)
    freezeColors
    
%Right Path Spheres
for i=1:sp_step:length(state_rgt)
    sphereE=ft2m((x_sp.*bubble) + state_rgt(i,1));
    sphereN=ft2m((y_sp.*bubble) + state_rgt(i,2));
    sphereU=ft2m((z_sp.*bubble) + state_rgt(i,3));
    
    [lat_sp, lon_sp, alt_sp] = enu2geodetic(sphereE,sphereN,sphereU,...
        lat_pos(c), lon_pos(c),h0,E,'radians');
    
    lat_sp=rad2deg(lat_sp);
    lon_sp=rad2deg(lon_sp);
    alt_sp = alt_sp - G_Height;
    h5=mesh(lon_sp,lat_sp,alt_sp);
    
end
    colormap(blue_sp)
    freezeColors
    
%Right-Up Path Spheres
for i=1:sp_step:length(state_rgt_up)
    sphereE=ft2m((x_sp.*bubble) + state_rgt_up(i,1));
    sphereN=ft2m((y_sp.*bubble) + state_rgt_up(i,2));
    sphereU=ft2m((z_sp.*bubble) + state_rgt_up(i,3));
    
    [lat_sp, lon_sp, alt_sp] = enu2geodetic(sphereE,sphereN,sphereU,...
        lat_pos(c), lon_pos(c),h0,E,'radians');
    
    lat_sp=rad2deg(lat_sp);
    lon_sp=rad2deg(lon_sp);
    alt_sp = alt_sp - G_Height;
    h6=mesh(lon_sp,lat_sp,alt_sp);
    
end
    colormap(purple_sp)
    freezeColors
    fprintf('\n'); 
    %Evaluates if a collision has occurred and breaks the loop if so
    if all(t_grnd_col_all(:) >= 0) 
        exitval = 1;
        break
    end
end
%exits the loop after all iterations are complete
exitval = 1;
break
end
xlabel('Longitude, deg','fontsize',14)
ylabel('Latitude, deg','fontsize',14)
zlabel('Altitude, m (MSL)','fontsize',14)
legend([h2 h3 h4 h5 h6],{'Forward Path','Left Path','Left-Up Path','Right Path','Right-Up Path'},'fontsize',14)
ylim([35.1 35.3])
xlim([-117.5 -117.33])

%This builds and plots the graphical path collision summary
    figure(2)
    %builds the spy matrix to be 1 for a collision 0 otherwise
    t_grnd_col_spy(t_grnd_col_spy(:) >= 0)  = 1;
    t_grnd_col_spy(t_grnd_col_spy(:) == -1) = 0;
    %flips columns for format [lft lft_up fwd rgt_up rgt]
    t_grnd_col_spy(:,[1 2 3 4 5])=t_grnd_col_spy(:,[2 4 1 5 3]); 
    iter_lim = (1:length(t_grnd_col_spy))';
    t_grnd_col_spy=[iter_lim,t_grnd_col_spy];
    t_grnd_col_spy = [t_grnd_col_spy(:,1),5*t_grnd_col_spy(:,2),4*t_grnd_col_spy(:,3),3*t_grnd_col_spy(:,4),2*t_grnd_col_spy(:,5),t_grnd_col_spy(:,6)];
    %Plots the spy matrix
    hg=plot(t_grnd_col_spy(:,1),t_grnd_col_spy(:,2:6),'rx','markersize',23,'linewidth',3);
    grid on
    set(gca,'xTick',1:1:length(t_grnd_col_spy))
    set(gca,'yTick',1:1:5)
    ylim([.5 5.5])
    set(gca,'yTickLabel',{'Rgt','Rgt-Up','Fwd','Lft-Up','Lft'})
    xlim([0 length(t_grnd_col_spy)+1])
    title('Graphical Path Collision Summary','fontsize',18,'fontweight','bold')
    xlabel('Iteration','fontsize',22,'fontweight','bold')
    ylabel('Path','fontsize',22)
    
    tilefigs