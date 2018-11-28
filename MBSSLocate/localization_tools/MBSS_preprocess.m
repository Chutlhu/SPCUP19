function [pairId, dMic, alpha, alphaSampled, tauGrid] = MBSS_preprocess(c, micPos, thetaGrid, phiGrid, alphaRes)

% Function MBSS_preprocess
% This function preprocess values needed to compute and aggregate
% the angular spectrum over all microphone pairs:
% - Combination of all microphone pair indexes
% - Angles corresponding to each direction {azimuth, elevation} for each microphone pair
% - Sampled angles array for each microphone pair
% - Corresponding TDOAs for each microphone pair
% This function depends only on the microphone positions and the search 
% angle grid boundaries & resolution (array of potential direction for sources).
%
% INPUT:
% c         : 1x1 , Speed of sound
% micPos    : N x 3 , cartesian coordinates of the N microphones
% thetaGrid : 1 x nDirection : Azimuth grid
% phiGrid   : 1 x nDirection : Elevation grid
% alphaRes  : 1x1, interpolation resolution
%
% OUTPUT:
% pairId       : nMicPair x 2, All microphone pair indexes
% dMic         : nMicPair x 1, distance between microphones for each pair
% alpha        : nMicPair x nDirection : Array of angles for each
%                microphone pair corresponding to all {theta, phi} to be
%                tested.
% alphaSampled : 1 x nMicPair cell array, each cell element contains the 
%                uniformly distributed angles to be tested for the 
%                corresponding pair
% tauGrid      : 1 x nMicPair cell array, each cell element contains the 
%                TDOA corresponding to the alphaSampled for each pair
%
% Version: v1.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2015 Ewen Camberlein and Romain Lebarbenchon
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)
% If you find it useful, please cite the following reference:
% http://bass-db.gforge.inria.fr/bss_locate/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of theta/phi combinations
nDirection = length(thetaGrid);

% Find all microphone pair indexes
nMic = size(micPos,1);
pairId = nchoosek(1:nMic,2);
nMicPair = size(pairId,1);

% Microphone direction vector (in xyz) for each pair
pfMn1n2 = micPos(pairId(:,1),:) - micPos(pairId(:,2),:);

% Microphone distance for each pair
dMic = sqrt(sum(pfMn1n2.^2,2));

% Convert all tuple {theta,phi} on the sphere grid in cartesian coordinates
Pjk = zeros(3,nDirection);
[Pjk(1,:), Pjk(2,:), Pjk(3,:)] = sph2cart(thetaGrid*pi/180, phiGrid*pi/180, 1);

% Note : Matlab specific method of computation to obtain alpha without loop - uses lots of RAM
% Duplicate all {theta,phi} coordinates for each pair of microphone 
Pjk_All = repmat(Pjk,[1 1 nMicPair]);
Pjk_All = permute(Pjk_All,[1 3 2]);

% Microphone direction vector duplicate for each {theta,phi} combination
Mn1n2_All = repmat(pfMn1n2',[1 1 nDirection]);

% alpha for one pair and one {theta,phi} is the angle formed between the microphone
% direction and the {theta,phi} direction - computation with dot product approach.
alpha = real(acosd(shiftdim(sum(Pjk_All.*Mn1n2_All),1)./repmat(dMic,[1 nDirection])));

% Compute 1D angles search grids and associated TDOA (Tau) search grids for each microphone pair  
% following search grid boundaries for each microphone pair is driven by
% the following fact : basic boundaries [0° 180°] for each pair could be
% adapted when the global search grid does not cover the entire space
% (leading to avoid useless curves computation and saving CPU time)
alphaSampled = cell(1,nMicPair);
tauGrid = cell(1,nMicPair);

for index = 1:nMicPair
    alphaSampled{index} = floor(min(alpha(index,:))/alphaRes) * alphaRes : alphaRes : ceil(max(alpha(index,:))/alphaRes) * alphaRes;
    tauGrid{index} = dMic(index)*cos(alphaSampled{index}.*pi/180)./c;
end

end
