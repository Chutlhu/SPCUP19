function s = MBSS_InputParam2Struct(angularSpectrumMeth,c,fftSize_sec,blockDurationUser_sec,blockOverlapUser_percent,pooling,azBound,elBound,gridRes,alphaRes,minAngle,nsrc,fs,normalizeSpecInst,specDisplay,enableWienerFiltering,wienerMode,freqRange,micPos,isArrayMoving,subArray,sceneTimeStamps,angularSpectrumDebug)
% Function MBSS_InputParam2Struct 
% Stack input parameters in a structure and apply default initialization if
% needed. In addition, compute the global search grid.
%
% Mandatory Inputs:
% nsrc                : 1 x 1, number of sources to be located
% fs                  : 1 x 1, sampling frequency (in hetz)
% micPos              : 3 x N x T, Cartesian positions of the N microphones over the time
% isArrayMoving       : 1 x 1, Flag to set if the microphone array is moving

% Optional Inputs (could be empty, in this case the default value will be applied):
%
% angularSpectrumMeth       : string, local angular spectrum method used :
%                            'GCC-PHAT'(default) 'GCC-NONLIN'  'MVDR' 'MVDRW' 'DS' 'DSW' 'DNM'
% c                         : 1 x 1, speed of sound (m/s). Default : 343 m/s
% wlen                      : 1 x 1, frame size / fft points. Default : 1024
% blockDurationUser_sec     : 1 x 1, Analysis window duration (in seconds). Default : entire signal. Note : rounded to have an integer number of frames in the block
% blockOverlapUser_percent  : 1 x 1, Overlap between successive blocks (in percent). Default 50%. Note : rounded to have an integer number of overlaped frames in the block
% pooling                   : string, pooling function (i.e. time integration method): 'max' (default) or 'sum'
% azBound                   : 1 x 2, azimuth boundaries (degrees). Default : [-180 180]
% elBound                   : 1 x 2, elevation boundaries (degrees). Default : [-90 90]
% gridRes                   : 1 x 1, sampling resolution applied to azimuth and elevation parameters (degrees); Default value: 1°
% alphaRes                  : 1 x 1, sampling resolution applied to each microphone pair referential (degrees); Default value: 5°
% minAngle                  : 1 x 1, minimum distance (degrees) between two peaks. Default value: 1°;
% normalizeSpecInst         : 1 x 1 , flag used to activate normalization of instantaneous local angular spectra. Default: 0
% specDisplay               : 1 x 1 , flag used to display the global angular spectrum and directions of sources found. Default: 0.
% enableWienerFiltering     : 1 x 1 , Flag used to apply Wiener filtering for signal or noise localization enhancement. Require an excerpt of noise signal. Default: 0.
% wienerMode                : string, enhancement mode if enableWienerFiltering set to 1. Could be 'Attenuation' or 'Emphasis'. Default: 'Attenuation'
% freqRange                 : 1 x 2 or [], vector containing the minimum and maximum frequency used by the aggregation function. If empty, all frequency bins are used. 
% subArray                  : 1 x K or [], Indices of the sub array microphones to use. [] to use all the microphones
% sceneTimeStamps           : T x 1 or [], Microphone position timeStamp, must be defined if T > 1
% angularSpectrumDebug      : 1 x 1 or [], Flag to enable additional plots to debug the angular spectrum aggregation. Default: disabled
%
% Output :
% s                   : struct of parameters containing the following fields
%   * nsrc                  : 1 x 1, number of sources to be located
%   * fs                    : 1 x 1, sampling frequency (in hetz)
%   * angularSpectrumMeth   : string, local angular spectrum method used :
%                             'GCC-PHAT' 'GCC-NONLIN'  'MVDR' 'MVDRW' 'DS' 'DSW' 'DNM'
%   * c                     : 1 x 1, speed of sound (m/s).
%   * useLinearTransform    : 1 x 1, flag, set to 1 for 'GCC-PHAT' or 'GCC-NONLIN' methods
%   * wlen                  : 1 x 1, frame size / fft points.
%   * blockProcessing       : true if block processing approach (i.e. blockDuration if neithor 'inf' nor 'empty'), false otherwise
%   * blockDuration_sec     : 1 x 1, Analysis window duration (in seconds). Default : entire signal. Note : rounded to have an integer number of frames in the block
%   * blockDuration_samples : 1 x 1, Analysis window duration (in seconds). Default : entire signal. Note : rounded to have an integer number of frames in the block
%   * blockDuration_frames  : 1 x 1, Analysis window duration (in frames). Default : entire signal. Note : rounded to have an integer number of frames in the block
%   * blockOverlap_percent  : 1 x 1, Overlap between successive blocks (in percent). Default 50%. Note : rounded to have an integer number of overlaped frames in the block
%   * blockStep_frames      : 1 x 1, block step used for successive blocks (in frames). Note : rounded to an integer number of frames
%   * blockStep_sec         : 1 x 1, block step used for successive blocks (in seconds). Note : rounded to an integer number of frames
%   * blockStep_samples     : 1 x 1, block step used for successive blocks (in samples). Note : rounded to an integer number of frames
%   * pooling               : string, pooling function (i.e. time integration method)
%   * alphaRes              : 1 x 1, sampling resolution applied to each microphone pair referential (degrees)
%   * minAngle              : 1 x 1, minimum distance (degrees) between two peaks.
%   * normalizeSpecInst     : 1 x 1 , flag used to activate normalization of instantaneous local angular spectra.
%   * specDisplay           : 1 x 1 , flag used to display the global angular spectrum and directions of sources found.
%   * azimuth               : nAzimuts x 1, vector of azimuths angles
%   * elevation             : 1 x nElevations, vector of elevation angles
%   * azimuthGrid            : 1 x nGrid, azimuth grid (nGrid = nAzimuts x nElevations)
%   * elevationGrid         : 1 x nGrid, elevation grid (nGrid = nAzimuts x nElevations)
%   * f                     : frequency axis (in Hertz)
%   * enableWienerFiltering  : 1 x 1, flag used to apply wiener filtering on input data. An example of noise must be provided.
%   * wienerMode : string, enhancement mode if enableWienerFiltering set to 1
%   * freqBins              : 1 x K containing the index of frequency bins used for the aggregation
%   * micPos                : 3 x K x T, Cartesian positions of the selected microphones over the time
%   * timeStamp             : T x 1 or [], Microphone position timeStamp, must be defined if T > 1
%   * isMoving              : 1 x 1, flag indicates if the array is moving or not
%   * arrayCentroid         : 3 x T, microphone array centroid over the time
%   * angularSpectrumDebug  : 1 x 1, Flag to enable additional plots to debug the angular spectrum aggregation
%
% Version : v2.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2018 Ewen Camberlein and Romain Lebarbenchon
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)
% If you find it useful, please cite the following reference:
% http://bass-db.gforge.inria.fr/bss_locate/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = struct;
%% Check mandatory params
if(isempty(nsrc))
    error('Number of estimated sources is a mandatory parameter, please put a value.');
else
    s.nsrc = nsrc;
end

if(isempty(fs))
    error('Sampling frequency is a mandatory parameter, please put a value.');
else
    s.fs = fs;
end

if(isempty(isArrayMoving))
    error('You must specify if the array is moving or not');
else
    s.isMoving = isArrayMoving;
end

if(isempty(micPos))
    error('You must specify the microphone positions');
else
    [dim1,N,T] = size(micPos);
    if(dim1~=3),error('The first dimension of micPos must be 3 - (x,y,z) coordinates');end
    if(isArrayMoving && (T == 1)), error('If the array is moving, the position must be defined for each time stamp'); end
    if(isArrayMoving && isempty(sceneTimeStamps)), error('If the array is moving, time stamps must be defined'); end
    if (T~=1 && (isempty(sceneTimeStamps) || length(sceneTimeStamps)~=T)), error('timeStamp must be correctly defined'); end
end

%% Check other params
% Check angularSpectrumMeth
if(isempty(angularSpectrumMeth))
    % default init : GCC-PHAT
    s.angularSpectrumMeth = 'GCC-PHAT';
elseif(~any(strcmp(angularSpectrumMeth, {'GCC-PHAT' 'GCC-NONLIN' 'MVDR' 'MVDRW' 'DS' 'DSW' 'DNM'})))
    error('Unknown local angular spectrum.');
else
    s.angularSpectrumMeth = angularSpectrumMeth;
end

% useLinearTransform  or not
if(any(strcmp(s.angularSpectrumMeth, {'GCC-PHAT','GCC-NONLIN'})))
    s.useLinearTransform = 1;
else
    s.useLinearTransform = 0;
end

if(isempty(c))
    s.c = 343;
else
    s.c = c;
end

% Init wlen (fftSize in samples)
if(isempty(fftSize_sec))
    fftSize_sec = 0.064; %default fftSize : 64 ms
end
% round wlen it to next pow 2 if necessary
requestedWlen = fftSize_sec * fs;
s.wlen = 2^nextpow2(requestedWlen);

% Check blockDuration
if(isempty(blockDurationUser_sec) || isinf(blockDurationUser_sec))
    % No block processing
    s.blockProcessing = false;
    s.blockDuration_frames  = []; % whole file entire frames
    s.blockDuration_sec     = []; % whole file entire frames
    s.blockDuration_samples = []; % whole file entire frames
else
    
    % test blockDurationUser conformity: 
    % should be a multiple of wlen/2 and handle at least 2 frames
    blockDurationUser_frames = blockDurationUser_sec*fs/s.wlen*2-1;
    blockDurationInternal_frames = max(2,floor(blockDurationUser_frames)); % At least 2 frames per block for loc estimation

    if(blockDurationInternal_frames == blockDurationUser_frames)
        % keep user choice
        blockDurationInternal_sec = blockDurationUser_sec;
    else
        % Modify blockDuration to have an integer number of frame
        blockDurationInternal_sec = (blockDurationInternal_frames + 1)*s.wlen/(2*fs);
        fprintf('[MBSS_InputParam2Struct] blockDuration has been modified to have an integer number of frame in the block. \n New blockDuration : %.3f sec\n',blockDurationInternal_sec);
    end
    
    s.blockProcessing = true;
    s.blockDuration_frames  = blockDurationInternal_frames;
    s.blockDuration_sec     = blockDurationInternal_sec;
    s.blockDuration_samples = blockDurationInternal_sec * fs;
end

% Check block overlap
if(s.blockProcessing)
    
    if(isempty(blockOverlapUser_percent))
        % Default value : No overlap
        s.blockOverlap_percent = 0;
        s.blockStep_frames     = s.blockDuration_frames;
        s.blockStep_sec        = s.blockDuration_sec;
        s.blockStep_samples    = s.blockDuration_samples;
    else
        
        % test blockStepUser conformity: should be a multiple of wlen/2
        blockStepUser_frames = s.blockDuration_sec*( 1 - blockOverlapUser_percent/100)*fs/s.wlen*2-1; 
        blockStepInternal_frames = max(0,floor(blockStepUser_frames));
        
        if(blockStepInternal_frames == blockStepUser_frames)
            % keep user choice
            blockOverlapInternal_percent = blockOverlapUser_percent;
        else
            % Round blockOverlap to overlap an integer number of frames
            blockOverlapInternal_percent = (1-(blockStepInternal_frames+1)*s.wlen/(2*fs*s.blockDuration_sec))*100;
            fprintf('[MBSS_InputParam2Struct] blockOverlap has been modified to have an integer number of step block frames. \n New blockOverlap : %.3f %%\n',blockOverlapInternal_percent);
        end
        
        s.blockOverlap_percent = blockOverlapInternal_percent;
        s.blockStep_frames     = blockStepInternal_frames;
        s.blockStep_sec        = (s.blockStep_frames + 1)*s.wlen/(2*fs);
        s.blockStep_samples    = s.blockStep_sec *fs;
    end
else
    % no block processing
    s.blockOverlap_percent =  0;
    s.blockStep_frames     =  [];
    s.blockStep_sec        =  [];
    s.blockStep_samples    =  [];
end

% Check pooling method
if(isempty(pooling))
    % Default pooling method : max
    s.pooling = 'max';
elseif(~any(strcmp(pooling, {'max' 'sum'})))
    error('Unknown pooling function.');
else
    s.pooling = pooling;
end

% Check alpha resolution
if(isempty(alphaRes))
    % Default alphaRes : 5°
    s.alphaRes = 5;
elseif(alphaRes < 0)
    error('Alpha resolution must be a positive value.');
else
    s.alphaRes = alphaRes;
end

% Check grid res
if(isempty(gridRes))
    gridRes = 1;
end

% Check min angle
if(isempty(minAngle))
    s.minAngle = 1;
elseif(minAngle < gridRes && nsrc>1)
    error('Minimum angle between two peaks has to be upper than arimut/elevation resolution');
else
    s.minAngle = minAngle;
end

%Check  normalizeSpecInst
if(isempty(normalizeSpecInst))
    s.normalizeSpecInst = 0;
elseif((nsrc == 1) && normalizeSpecInst)
    warning('Use of instantaneous local angular spectra normalization with one source to be located is unnecessary. Switch off normalization.');
    s.normalizeSpecInst = 0;
else
    s.normalizeSpecInst = normalizeSpecInst;
end

% check specDisplay
if(isempty(specDisplay))
    s.specDisplay = 0;
else
    s.specDisplay = specDisplay;
end

% check enableWienerFiltering
if(isempty(enableWienerFiltering) || enableWienerFiltering == 0)
    s.enableWienerFiltering = 0;
else
    if(s.useLinearTransform == 0)
        error('Could not apply wiener filtering if the method does not use a linear transform. Change your method to one of the following: GCC-PHAT or GCC-NONLIN');
    else
        s.enableWienerFiltering = enableWienerFiltering;
    end
end

% check wienerMode
if(isempty(enableWienerFiltering))
    s.wienerMode = ''; % not used
else
    if(isempty(wienerMode))
        s.wienerMode = 'Attenuation'; % default value
    else
        if(~any(strcmp(wienerMode, {'Attenuation' 'Emphasis'})))
            error('Unknown wienerMode');
        else
            s.wienerMode = wienerMode;
        end
    end
end

% angularSpectrumDebug
if(isempty(angularSpectrumDebug))
    s.angularSpectrumDebug = 0;
else
    s.angularSpectrumDebug = angularSpectrumDebug;
end
    

% Compute the frequency axis
s.f = s.fs/s.wlen*(1:s.wlen/2).';

% frequency bins used for the aggregative part
if(isempty(freqRange))
    % All the bins are used
    s.freqBins = 1:length(s.f);
elseif(freqRange(1) < 0 || freqRange(2) > s.fs/2)
    error('Frequency range must be between 0Hz and Fs/2');
else
    binMin = find(s.f >= freqRange(1),1,'first');
    binMax = find(s.f<freqRange(2),1,'last');
    s.freqBins = binMin:binMax;
end
%% Compute the search grid
if(isempty(azBound))
    azBound = [-180 180];
elseif((length(azBound) == 1) && azBound >= -90 && azBound <= 90)
    azBound = [azBound,azBound];
elseif(length(azBound) == 2 && azBound(1) >= -180 && azBound(2) <= 180 && azBound(1)<=azBound(2))
    % nothing to do
else
    error('Azimut boundaries are bad filled. Azimut boundaries could be:\n - One scalar value to locate at a specific azimuth\n - A vector of two ascending values between -/+ 180°');   
end

if(isempty(elBound))
    elBound = [-90 90];
elseif(length(elBound) == 1 && elBound >= -90 && elBound <= 90)
    elBound = [elBound,elBound];
elseif(length(elBound) == 2 && elBound(1) >= -90 && elBound(2) <= 90 && elBound(1)<=elBound(2))
    % nothing to do
else
    error('Elevation boundaries are bad filled. Elevation boundaries could be:\n - One scalar value to locate at a specific elevation\n - A vector of two ascending values between -/+ 90°');   
end

if(length(unique(elBound)) == 1 && length(unique(azBound)) == 1)
    error('You can not fixed the boundaries for azBound and elBound to an unique value');
end

s.azimuth = (azBound(1) : gridRes : azBound(2))';
s.elevation   = (elBound(1) : gridRes : elBound(2));
nAz = length(s.azimuth);
nEl = length(s.elevation);
s.azimuthGrid = repmat(s.azimuth,nEl,1)';
s.elevationGrid   = (reshape(repmat(s.elevation,nAz,1),1,nAz*nEl));

if (~isempty(subArray))
    if( (length(subArray) == 1 || length(subArray) > N) || min(subArray) <= 0 || max(subArray) > N)
        error('subArray is badly defined');        
    end
else
    subArray = 1:N; % all the mics
end

% Fill the structure
s.micPos = micPos(:,subArray,:);
s.timeStamp = sceneTimeStamps;
s.isMoving = isArrayMoving;

% Compute the microphone array centroid
s.arrayCentroid = squeeze(mean(s.micPos,2));
end