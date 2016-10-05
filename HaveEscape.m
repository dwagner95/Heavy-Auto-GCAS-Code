function [active,pathCollision,preTimer,pathTimers,loopTimer,postTimer,fullTimer] = HaveEscape(...
    latitude, longitude, altitude, course, gamma ...
    ,odeTimeVector,odeCmdVectorTime,odeCmdMatrixNz,odeCmdMatrixPhi ...
    ,dtedLatMatrix,dtedLonMatrix,dtedHotMatrix ...
    ,pathOnly)
%#codegen

persistent pActivePath;
persistent pPathCollisionVector;

if isempty(pActivePath) 
    pActivePath = 0;
end

if isempty(pPathCollisionVector)
    pPathCollisionVector=zeros(5,1);
end

if ( pathOnly > 0 && pathOnly <=5)
    pActivePath = pathOnly;
end

if ( pActivePath > 0 && pActivePath <=5 )
    active = pActivePath;
    fullTimer = 0;
    preTimer = 0;
    pathTimers = zeros(5,1);
    loopTimer = 0;
    postTimer = 0;
    pathCollision = pPathCollisionVector;
    return;
end

coder.extrinsic('tic');
coder.extrinsic('toc');
fullTic = uint64(zeros);
fullTic = tic;




V=523.22;                       %Constant Velocity, ft/s (310 kts)
altitudeM = ft2m(altitude);               %MSL Altitude as read from Altimeter (ft)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%No User Inputs Below This Line%%%%%%%%%%%%%%%%%%%%%%%%%%

global gGHeight;

altitudeM = altitudeM + gGHeight;           %Current WGS-84 HAE A/C Altitude (m)

dtedLatVector = dtedLatMatrix(:,1);
dtedLonVector = dtedLonMatrix(1,:)';

hot0 = interp2(dtedLatVector,dtedLonVector,dtedHotMatrix',...
    latitude,longitude);

% if ~isfinite(hot0)
%     hot0 = altitudeM;
% end

[xM,yM,zM] = my_geodetic2enu(...
    latitude,longitude,altitudeM,...
    latitude,longitude,hot0);


% Initialize A/C state values for Propagation with ODE45
xFt = m2ft(xM);     % Initial x position, ft
yFt = m2ft(yM);     % Initial y position, ft
zFt = m2ft(zM);     % Initial z position, ft (Altitude)
%chiDeg = compass2pol(headingAngleDeg);    % Initial heading angle
coder.extrinsic('wrapTo360');
chiDeg = zeros;
chiDeg = wrapTo360(90-course);

%Initial conditions column vector [x;y;z;gam;chi] [ft,ft,ft,deg,deg]
initialConditions= [xFt;yFt;zFt;deg2rad(gamma);deg2rad(chiDeg)];

global gNumberOfPaths;
global gFeetPerDegLonPoly;

if gNumberOfPaths > 5
    gNumberOfPaths = 5;
end

% Allocate the result
% numberOfPaths = gNumberOfPaths;
% assert(odeCmdTimeLength < 1000);
% assert(numberOfPaths<6);
% aircraftProjectedPositionM = zeros(...
%     odeCmdTimeLength,...
%     3,...
%     numberOfPaths...
%     );

% Allocate result
pathCollisionVector = zeros(5,1);
pathTimerTemp = zeros(5,1);

global gBubble;
global gM;
global gPredictionTimeBuffer;

preTimer = zeros;
preTimer = toc(fullTic);

loopTic = uint64(zeros);
loopTic = tic;

for i=1:gNumberOfPaths
    
    pathTic = uint64(zeros);
    pathTic = tic;

    
    [state] = myOde( ...
        odeTimeVector,...
        initialConditions,...
        gPredictionTimeBuffer,...
        [odeCmdVectorTime odeCmdMatrixNz(:,i)],...
        [odeCmdVectorTime -odeCmdMatrixPhi(:,i)],...
        V );
    
    aircraftProjectedPositionM = ft2m(state(:,1:3));
    
% end
% 
% for i=1:gNumberOfPaths
    
    %Convert propagated position back into Lat and Lon
    [pathVectorLatDeg,pathVectorLonDeg,pathVectorAlt]=my_enu2geodetic(...
        aircraftProjectedPositionM(:,1),aircraftProjectedPositionM(:,2),aircraftProjectedPositionM(:,3),...
        latitude, longitude,hot0);
    %pathLatVectorRad=rad2deg(pathLatVectorDeg);
    %pathLonVectorRad=rad2deg(pathLonVectorDeg);
    
    
    
    
    % Logic to check the solution for terrain collision
    minDistanceBetweenAircraftAndTerrainM = zeros(1,length(odeTimeVector));
    for j=1:length(odeTimeVector)
        
        %Finds the distance between a deg of longitude at a specific Lat in km
        kmPerDegLon = ft2m(polyval(gFeetPerDegLonPoly,pathVectorLatDeg(j)))/1000;
        
        %Conversion from km to degrees for longitude distance
        bubbleRadiusinDegLon = inv(kmPerDegLon)*inv(3280.84)*gBubble; %deg Long at Lat position
        bubbleRadiusinDegLat = km2deg(gBubble);
        
        %finds the column indices Longitude DTED Posts
        dtedIndicesWithinBubbleLon=find( dtedLonVector <= (pathVectorLonDeg(j)+bubbleRadiusinDegLon)...
                                       & dtedLonVector >= (pathVectorLonDeg(j)-bubbleRadiusinDegLon) );
        
        %Finds the row indices Latitude DTED Posts
        dtedIndicesWithinBubbleLat=find( dtedLatVector <= (pathVectorLatDeg(j)+bubbleRadiusinDegLat)...
                                       & dtedLatVector >= (pathVectorLatDeg(j)-bubbleRadiusinDegLat) );
        
        %Breaks current loop iteration if no DTED Posts found
        if isempty(dtedIndicesWithinBubbleLon)==1 || isempty(dtedIndicesWithinBubbleLat)==1
            minDistanceBetweenAircraftAndTerrainM(j) = 8888;
            continue
        end
        
        
        %Finds Long/Lat values from complete Long/Lat Matrices
        localDtedLonMatrix = dtedLonMatrix(dtedIndicesWithinBubbleLat,dtedIndicesWithinBubbleLon);
        localDtedLatMatrix = dtedLatMatrix(dtedIndicesWithinBubbleLat,dtedIndicesWithinBubbleLon);
        localDtedHotMatrix = dtedHotMatrix(dtedIndicesWithinBubbleLat,dtedIndicesWithinBubbleLon); %(long,lat) because Zgrid is transposed
        
        %Build lat_lon_alt DTED posts for this timestep
        kount = numel(localDtedLonMatrix);
        assert(kount<722);
        localDtedPostList = zeros(kount,3);            
        for k=1:kount
            localDtedPostList(k,1:3)=[localDtedLatMatrix(k),localDtedLonMatrix(k),localDtedHotMatrix(k)];
        end
        
        %convert DTED to ENU
        [x,y,z]=my_geodetic2enu(...
            localDtedPostList(:,1),localDtedPostList(:,2),localDtedPostList(:,3),...
            latitude, longitude, hot0);
            %Answer in meters
        
        %Build matrices required to check if DTED are inside bubble
        dtedPostPositionMatrixM = [x,y,z];
        
        %Checks if DTED Post height is inside sphere XMX'
        xyzDtedPostBubbleOverlayLength = length(x);
        assert(xyzDtedPostBubbleOverlayLength<100);
        aircraftDistanceToDtedPost = zeros(xyzDtedPostBubbleOverlayLength,1);
        for k=1:xyzDtedPostBubbleOverlayLength
            aircraftDistanceToDtedPost(k) = ...
                (dtedPostPositionMatrixM(k,:)- aircraftProjectedPositionM(j,:))...
               *gM...
               *(dtedPostPositionMatrixM(k,:)- aircraftProjectedPositionM(j,:))';
        end
        
        %Finds the smallest distance for a given timestep
        minDistanceBetweenAircraftAndTerrainM(j) = min(aircraftDistanceToDtedPost);
    end
    
    indexOfGroundCollisions = find(minDistanceBetweenAircraftAndTerrainM <=1);
    
    %Calculates the fwd time of a collision occurring.
    
    if isempty(indexOfGroundCollisions) == 0 %Meaning: not empty
        timeOfFirstGroundCollision = odeTimeVector(indexOfGroundCollisions(1));
        pathCollisionVector(i) = timeOfFirstGroundCollision(1);
    else
        %Fills an empty matrix with -1 to delineate that no collision occured
        pathCollisionVector(i) = -1;
    end
    
    pathTimerTemp(i) = toc(pathTic);
    
end

loopTimer = zeros;
loopTimer = toc(loopTic);



    postTic = uint64(zeros);
    postTic = tic;
    
% 
% pathCollision = pPathCollisionVector;
% 


% All paths collide
if isempty(find(pathCollisionVector<0,1))
    % Look at the last collision vector
%     clearPaths = find(pPathCollisionVector<0);
%     pActivePath = clearPaths(1);
    
    % OR
    
    % Pick the largest
    largest = find(pathCollisionVector==max(pathCollisionVector));
    pActivePath = largest(1);
end


% switch length(find(pathCollisionVector<0))
%     case 0 % If all paths collide
%         pActivePath = 3; % Pick LARGEST
%     case 1 % If one escape route
%         temp = find(pathCollisionVector<0);
%         pActivePath = temp(1);
% end

pPathCollisionVector =   pathCollisionVector;
pathCollision = pPathCollisionVector;


% if all(pathCollisionVector(:) >= 0)
%     if t_grnd_col_fwd >= t_grnd_col_lft && t_grnd_col_fwd >= t_grnd_col_rgt...
%             && t_grnd_col_fwd >= t_grnd_col_lft_up && t_grnd_col_fwd >= t_grnd_col_rgt_up
%         pActivePath = 1;
%     elseif t_grnd_col_lft > t_grnd_col_fwd && t_grnd_col_lft > t_grnd_col_rgt...
%             && t_grnd_col_lft > t_grnd_col_lft_up && t_grnd_col_lft > t_grnd_col_rgt_up
%         pActivePath = 2;
%     elseif t_grnd_col_rgt > t_grnd_col_fwd && t_grnd_col_rgt > t_grnd_col_lft...
%             && t_grnd_col_rgt > t_grnd_col_lft_up && t_grnd_col_rgt > t_grnd_col_rgt_up
%         pActivePath = 3;
%     elseif t_grnd_col_lft_up > t_grnd_col_fwd && t_grnd_col_lft_up >= t_grnd_col_rgt...
%             && t_grnd_col_lft_up >= t_grnd_col_lft && t_grnd_col_lft_up >= t_grnd_col_rgt_up
%         pActivePath = 4;
%     elseif t_grnd_col_rgt_up > t_grnd_col_fwd && t_grnd_col_rgt_up >= t_grnd_col_rgt...
%             && t_grnd_col_rgt_up >= t_grnd_col_lft && t_grnd_col_rgt_up > t_grnd_col_lft_up
%         pActivePath = 5;
%     end
% end
% 

active = pActivePath;

fullTimer = zeros;
fullTimer = toc(fullTic);

postTimer = zeros;
postTimer = toc(postTic);



pathTimers = pathTimerTemp;