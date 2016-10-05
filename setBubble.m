function [ M ] = setBubble( bubbleFt )
%SETBUBBLE Summary of this function goes here
%   Detailed explanation goes here

%%%%%%%%%%%%%%%%%%%%%Build Terrain Map From DTED%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Map Read: GCAS Mtn in R-2508 at Edwards AFB CA

% latlim = [35.1 35.3];
% lonlim = [-117.5 -117.3];
% 
% % if bubbleFt < 300
% %     [Z, refvec] = dted('n35.dt2', 1, latlim, lonlim); % Elevation (m)
% % else
%     [Z, refvec] = dted('n35.dt1', 1, latlim, lonlim); % Elevation (m)
% % end
% 
% cellsperdegree = refvec(1);
% step = 1/cellsperdegree;
% 
% Zlength = length(Z);            %meters
% latN = refvec(2);     %pulls the northern lat position from refvec
% lonW = refvec(3);     %pulls the Western lon position from refvec
% lonE = (Zlength-1)*step+lonW;  %calculates the opposite corner
% latS = (Zlength-1)*(-step)+latN; %calculates the opposite corner
% 
% xVector = lonW:step:lonE;
% xLongMatrix = repmat(xVector,Zlength,1); %deg
% 
% yVector = latS:step:latN;
% yLatMatrix = repmat(yVector,Zlength,1)'; %deg
% 
% zGrid = (Z + gHeight)'; %Turns all DTED altitude values into HAE
% %haeGroundHeightInterpolant = griddedInterpolant(xGrid,yGrid,zGrid,'spline'); %HAE Altitudes
% 

bubble_meters = ft2m(bubbleFt);%Puts bubble in meters for consistency
R2 = bubble_meters^2;          %Squared once for speed (m^2)
M = diag(1./[R2 R2 R2]);       %For transformation matrix (m^2)


end

