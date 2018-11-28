function [azEst, elEst, blockTimestamps,sElapsTime,figHandle] = MBSS_locate_spec(x,wienerRefSignal,sMBSSParam)
%% Function MBSS_locate_spec
% Estimate source localizations in a multichannel convolutive mixture using
% an angular spectrum based approach
%
% Inputs:
% x             : nsampl x nchan, matrix containing nchan time-domain mixture
%                 signals with nsampl samples
% wienerRefSignal : nsampl x nchan, matrix containing a nchan time-domain
%                 signal to be attenuated or emphazed into the mixture x.
%                 Only used by enabling the Wiener filtering parameter.
%
% sMBSSParam    : structure, MBSS parameters structure containing the following fields
%   * nsrc                  : 1 x 1, number of sources to be located
%   * fs                    : 1 x 1, sampling frequency (in hetz)
%   * angularSpectrumMeth   : string, local angular spectrum method used :
%                             'GCC-PHAT' 'GCC-NONLIN'  'MVDR' 'MVDRW' 'DS' 'DSW' 'DNM'
%   * c                     : 1 x 1, speed of sound (m/s).
%   * useLinearTransform    : 1 x 1, flag, set to 1 for 'GCC-PHAT' and 'GCC-NONLIN'
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
%   * azimuthGrid           : 1 x nGrid, azimuth grid (nGrid = nAzimuts x nElevations)
%   * elevationGrid         : 1 x nGrid, elevation grid (nGrid = nAzimuts x nElevations)
%   * f                     : frequency axis (in Hertz)
%   * enableWienerFiltering : 1 x 1, flag used to apply wiener filtering on input data. An example of noise must be provided.
%   * wienerMode            : string, enhancement mode if enableWienerFiltering set to 1
%   * freqBins              : 1 x K containing the index of frequency bins used for the aggregation
%   * micPos                : 3 x K x T, Cartesian positions of the selected microphones over the time
%   * timeStamp             : T x 1 or [], Microphone position timeStamp, must be defined if T > 1
%   * isMoving              : 1 x 1, flag indicates if the array is moving or not
%   * arrayCentroid         : 3 x T, microphone array centroid over the time
%   * angularSpectrumDebug  : 1 x 1, Flag to enable additional plots to debug the angular spectrum aggregation
%
% Outputs:
% azEst           : nBlocks x nsrc, vector of estimated azimuths (in degrees) on each block. NaN value if not found
% elEst           : nBlocks x nsrc, vector of estimated elevations (in degrees) on each block. NaN value if not found
% blockTimestamps : nBlocks x 1, Time stamp of the block of frames (middle of the block)
% sElapsTime      : structure, This structure gives the following processing time estimates (in seconds) :
%   * xTFComp  : 1 x 1, processing time of the Time-Frequency transform on input x
%   * wienerRefSignalTFComp : 1 x 1, processing time of the Time-Frequency transform on input wienerRefSignal
%   * wienerRefSignalCovComp: 1 x 1, processing time of the noise covariance matrix computation
%   * blockLoc : nBlocks x 1, processing time of the localization on each block
% figHandle        : 1 x 1, displayed figure handle  
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


%% Check the inputs
[nSamples,nChan]=size(x);
[~,nMic,T] = size(sMBSSParam.micPos);
if nChan>nSamples, error('The input signal must be in columns.'); end
if nChan~=nMic, error('Number of microphones and number of signal channels must be the same'); end

fprintf('Input signal duration: %.02f seconds\n',nSamples/sMBSSParam.fs);

if(~sMBSSParam.isMoving)
    % Static microphone array => MBSS_preprocess is called once
    MBSS_PREPROCESS_ONCE = true;
else
    % Moving microphone array => MBSS_preprocess is called for each chunk
    % of frames
    MBSS_PREPROCESS_ONCE = false;
end

% Init the sElapsTime output structure
sElapsTime =  struct;
tstart = tic;
%% Compute the Time-Frequency Representation
if(sMBSSParam.useLinearTransform)
    % Linear transform
    [X,startSample,~] = MBSS_stft_multi(x.',sMBSSParam.wlen);
    X = X(2:end,:,:);  
else
    % Quadratic transform
    [X,startSample,~] = MBSS_qstft_multi(x,sMBSSParam.fs,sMBSSParam.wlen,8,2);
    X = permute(X(:,:,2:end,:),[3 4 1 2]);
end
sElapsTime.xTFComp = toc(tstart);

nframe = size(X,2);

%% Deal the block processing approach
if(sMBSSParam.blockProcessing)
    % Only keep full block of frames
    frameStart = 1 : sMBSSParam.blockStep_frames : nframe;
    frameEnd = sMBSSParam.blockDuration_frames : sMBSSParam.blockStep_frames : nframe;
    nblocks = min(length(frameStart),length(frameEnd));
    frameStart = frameStart(1:nblocks);
    frameEnd = frameEnd(1:nblocks);
    blockTimestamps = ((startSample(frameEnd) + startSample(frameStart))/2)./sMBSSParam.fs;
else
    frameStart = 1;
    frameEnd = nframe;
    nblocks = 1;
    blockTimestamps = ((startSample(frameEnd) + startSample(frameStart))/2)./sMBSSParam.fs;
end

%% Prepare outputs
% Output matrix init
azEst = nan(nblocks,sMBSSParam.nsrc); % estimated azimuth in degrees
elEst = nan(nblocks,sMBSSParam.nsrc); % estimated elevation in degrees

%% Prepare displaying
if(sMBSSParam.specDisplay)
    figHandle = figure;
else
    figHandle = -1;
end

sElapsTime.blockLoc = nan(nblocks,1);% elapsed time in seconds

%% Compute the noise covariance matrix if needed
sElapsTime.wienerRefSignalTFComp = 0;
sElapsTime.wienerRefSignalCovComp = 0;
if(sMBSSParam.enableWienerFiltering == 1)
    % wiener filtering
   
    if(isempty(wienerRefSignal))
        error('In order to apply wiener filtering on mixture, you must provide wienerRefSignal');
    else
        % compute the STFT on noise input samples
        tstart = tic;
        Xn = MBSS_stft_multi(wienerRefSignal',sMBSSParam.wlen);
        Xn = Xn(2:end,:,:);
        sElapsTime.wienerRefSignalTFComp = toc(tstart);
        
        % Compute the covariance on noise
        tstart = tic;
        Rnn = MBSS_covariance(Xn);
        sElapsTime.wienerRefSignalCovComp = toc(tstart);
        
    end

    fprintf('Elaps. time of noise covariance matrix computation : %f sec.\n',sElapsTime.wienerRefSignalCovComp + sElapsTime.wienerRefSignalTFComp);
else
    Rnn = [];
end
   
%% Process each block of frames
sElapsTime.blockLoc = nan(nblocks,1);

for block_idx = 1 : nblocks
    tstart = tic;
    fprintf('Block %d / %d \n',block_idx,nblocks);
    if(~MBSS_PREPROCESS_ONCE || (block_idx == 1))
        
        if(T > 1)
            % Find the nearest microphone position :
            diff = abs(bsxfun(@minus,blockTimestamps(block_idx), sMBSSParam.timeStamp));
            [~,closest_idx] = min(diff);
            
            % Get the microphone position at the closest temporal index
            micPos = shiftdim(sMBSSParam.micPos(:,:,closest_idx));

        else
            micPos = sMBSSParam.micPos;
        end
        
        % Create the aggregationParam structure
        aggregationParam = struct;
        [aggregationParam.pairId, aggregationParam.d, aggregationParam.alpha, aggregationParam.alphaSampled, aggregationParam.tauGrid] = ...
            MBSS_preprocess(sMBSSParam.c, micPos', sMBSSParam.azimuthGrid, sMBSSParam.elevationGrid, sMBSSParam.alphaRes);
        
        aggregationParam.nGrid = length(sMBSSParam.azimuthGrid);
        aggregationParam.nPairs = size(aggregationParam.pairId,1);
        aggregationParam.c = sMBSSParam.c;
        aggregationParam.fs = sMBSSParam.fs;
        aggregationParam.azimuth = sMBSSParam.azimuth;
        aggregationParam.elevation = sMBSSParam.elevation;
    end
    
    % Apply wiener filtering on X if requested
    if(sMBSSParam.enableWienerFiltering)
        X_current = MBSS_wiener(X(:,frameStart(block_idx):frameEnd(block_idx),:,:,:),Rnn,sMBSSParam.wienerMode);
    else
        X_current = X(:,frameStart(block_idx):frameEnd(block_idx),:,:,:);
    end
    
    % Compute the angular spectrum
    specInst = MBSS_computeAngularSpectrum(sMBSSParam.angularSpectrumMeth,sMBSSParam.angularSpectrumDebug,aggregationParam,X_current,sMBSSParam.f,sMBSSParam.freqBins);
    
    % Normalize instantaneous local angular spectra if requested
    if(sMBSSParam.normalizeSpecInst)
        [~,nFrames,~] = size(specInst);
        for i=1:nFrames
            minVal = min(min(specInst(:,i)));
            specInst(:,i)=(specInst(:,i) - minVal)/ max(max(specInst(:,i)- minVal));
        end
    end
    
    % Pooling
    switch sMBSSParam.pooling
        case 'max'
            specGlobal = shiftdim(max(specInst,[],2));
        case 'sum'
            specGlobal = shiftdim(sum(specInst,2));
    end
    
    [pfEstAngles,figHandle] = MBSS_findPeaks2D(figHandle,specGlobal, sMBSSParam.azimuth, sMBSSParam.elevation, sMBSSParam.azimuthGrid, sMBSSParam.elevationGrid, sMBSSParam.nsrc, sMBSSParam.minAngle, sMBSSParam.angularSpectrumMeth, sMBSSParam.specDisplay);
    nSrcFound = size(pfEstAngles,1);

    azEst(block_idx,1:nSrcFound) = pfEstAngles(:,1)';
    elEst(block_idx,1:nSrcFound) = pfEstAngles(:,2)';
    
    fprintf('Src: %s\n',num2str(1:nSrcFound));
    fprintf('Az.: %s\n',num2str(round(azEst(block_idx,1:nSrcFound))));
    fprintf('El.: %s\n\n',num2str(round(elEst(block_idx,1:nSrcFound))));
    
    sElapsTime.blockLoc(block_idx) = toc(tstart);
    fprintf('Elaps. time : %f sec.\n',sElapsTime.blockLoc(block_idx));
    
end
end

%% PEAKS FINDING METHODS
function [pfEstAngles,figHandle] = MBSS_findPeaks2D(figHandle,ppfSpec, piAzimuts, piElevations, piAzimutGrid, piElevationGrid, iNbSources, iMinAngle,angularSpectrumMeth, bDisplayResults)

% Function MBSS_findPeaks2D
%
% This function search peaks in computed angular spectrum with respect to iNbSources and iMinAngle.
%
% INPUT:
% ppfSpec          : 1 x iNGridDirection : 1D angular spectrum
% piAzimuts        : 1 x iNbThetas : Azimuth sampled values
% piElevations     : 1 x iNbPhis : Elevation sampled values
% piAzimutGrid     : 1 x iNGridDirection : Azimuth grid
% piElevationGrid  : 1 x iNGridDirection : Elevation grid
% iNbSources       : 1 x 1 : Number of sources to be found
% iMinAngle        : 1 x 1 : Minimum angle between two sources
% bDisplayResults  : 1 x 1 : display global 2D angular spectrum with sources found
% figHandle        : 1 x 1 
% OUTPUT:
% pfEstAngles : iNbSourcesFound x 2 : azimuth and elevation for all sources found
%
% Version: v2.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2018 Ewen Camberlein and Romain Lebarbenchon
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)
% If you find it useful, please cite the following reference:
% http://bass-db.gforge.inria.fr/bss_locate/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iNbThetas = length(piAzimuts);
iNbPhis = length(piElevations);

if(iNbThetas == 1 || iNbPhis == 1)
    % 1D search
    [~, ind] = findpeaks(ppfSpec, 'minpeakdistance',iMinAngle, 'sortstr','descend');
    iNbSourcesFound = min(iNbSources,length(ind));
    
    if(iNbThetas == 1)
        pfEstAngles = [piAzimuts.*ones(iNbSourcesFound,1) piElevationGrid(ind(1:iNbSourcesFound))'];
    else
        pfEstAngles = [piAzimutGrid(ind(1:iNbSourcesFound))' piElevations.*ones(iNbSourcesFound,1)];
    end
    
    % Display
    if (bDisplayResults)
        figure(figHandle);
        if(iNbThetas == 1)
            plot(piElevationGrid,ppfSpec);
            hold on
            %display sources found
            plot(piElevationGrid(ind(1:iNbSourcesFound)),ppfSpec(ind(1:iNbSourcesFound)),'*k','MarkerSize',15,'linewidth',1.5);
            xlabel('\phi (degrees)');
            ylabel('Angular spectrum');
            hold off;
            title(['\Phi^{' angularSpectrumMeth '}(\theta,\phi) for \theta = ' num2str(piAzimuts) '  |  markers : sources found ']);
        else
            plot(piAzimutGrid,ppfSpec);
            hold on
            %display sources found
            plot(piAzimutGrid(ind(1:iNbSourcesFound)),ppfSpec(ind(1:iNbSourcesFound)),'*k','MarkerSize',15,'linewidth',1.5);
            xlabel('\theta (degrees)');
            ylabel('Angular spectrum');
            hold off;
            title(['\Phi^{' angularSpectrumMeth '}(\theta,\phi) for \phi = ' num2str(piElevations) '  |  markers : sources found ']);
        end
    end
    % 
else
% Convert angular spectrum in 2D
ppfSpec2D = (reshape(ppfSpec,iNbThetas,iNbPhis))';


% search all local maxima (local maximum : value higher than all neighborhood values)
% some alternative implementations using matlab image processing toolbox are explained here :
% http://stackoverflow.com/questions/22218037/how-to-find-local-maxima-in-image)

% Current implementation uses no specific toolbox. Explanations can be found with following link :
% http://stackoverflow.com/questions/5042594/comparing-matrix-element-with-its-neighbours-without-using-loop-in-matlab
% All values of flat peaks are detected as peaks with this implementation :
ppfPadpeakFilter = ones(size(ppfSpec2D,1)+2,size(ppfSpec2D,2)+2) * -Inf;
ppfPadpeakFilter(2:end-1,2:end-1) = ppfSpec2D;

% Find peaks : compare values with their neighbours
ppiPeaks = ppfPadpeakFilter(2:end-1,2:end-1) >= ppfPadpeakFilter(1:end-2,2:end-1) & ... % top
    ppfPadpeakFilter(2:end-1,2:end-1) >= ppfPadpeakFilter(3:end,  2:end-1) & ... % bottom
    ppfPadpeakFilter(2:end-1,2:end-1) >= ppfPadpeakFilter(2:end-1,1:end-2) & ... % right
    ppfPadpeakFilter(2:end-1,2:end-1) >= ppfPadpeakFilter(2:end-1,3:end)   & ... % left
    ppfPadpeakFilter(2:end-1,2:end-1) >= ppfPadpeakFilter(1:end-2,1:end-2) & ... % top/left
    ppfPadpeakFilter(2:end-1,2:end-1) >= ppfPadpeakFilter(1:end-2,3:end)   & ... % top/right
    ppfPadpeakFilter(2:end-1,2:end-1) >= ppfPadpeakFilter(3:end,  1:end-2) & ... % bottom/left
    ppfPadpeakFilter(2:end-1,2:end-1) >= ppfPadpeakFilter(3:end,  3:end);        % bottom/right

% number of local maxima
iNbLocalmaxima = sum(sum(ppiPeaks));

% local maxima with corrresponding values
ppfSpec2D_peaks = (ppfSpec2D - min(min(ppfSpec2D))) .* ppiPeaks; % substract min value : avoid issues (when sorting peaks) if some peaks values are negatives

% sort values of local maxima
pfSpec1D_peaks= reshape(ppfSpec2D_peaks',1,iNbPhis*iNbThetas);
[~,piIndexPeaks1D] = sort(pfSpec1D_peaks,'descend');

piEstSourcesIndex = piIndexPeaks1D(1);  % first source is the global maximum (first one in piSortedPeaksIndex1D)
index = 2; % search index in piSortedPeaksIndex1D
iNbSourcesFound = 1; % set to one as global maximum is already selected as source

%Filter the list of peaks found with respect to minAngle parameter
while (iNbSourcesFound < iNbSources && index <= iNbLocalmaxima)
    
    bAngleAllowed = 1;
    % verify that current direction is allowed with respect to minAngle and sources already selected
    for i = 1:length(piEstSourcesIndex)
        
        % distance calculated using curvilinear abscissa (degrees) - ref. : http://geodesie.ign.fr/contenu/fichiers/Distance_longitude_latitude.pdf
        dist=acosd(sind(piElevationGrid(piEstSourcesIndex(i)))*sind(piElevationGrid(piIndexPeaks1D(index)))+cosd(piElevationGrid(piEstSourcesIndex(i)))*cosd(piElevationGrid(piIndexPeaks1D(index)))*cosd(piAzimutGrid(piIndexPeaks1D(index))-piAzimutGrid(piEstSourcesIndex(i))) );
        
        if(dist <iMinAngle)
            bAngleAllowed =0;
            break;
        end
    end
    
    % store new source
    if(bAngleAllowed)
        piEstSourcesIndex = [piEstSourcesIndex,piIndexPeaks1D(index)];
        iNbSourcesFound = iNbSourcesFound +1;
    end
    
    index = index + 1;
end

pfEstAngles = [piAzimutGrid(piEstSourcesIndex)' piElevationGrid(piEstSourcesIndex)'];

%% Display results
if (bDisplayResults)
    figure(figHandle);
    colormap(jet); % bleu jaune rouge
    imagesc(piAzimuts,piElevations,ppfSpec2D);
    set(gca,'YDir','normal');
    hold on;
    
    %display sources found
    for i =1:length(pfEstAngles(:,1))
        handle=plot(pfEstAngles(i,1),pfEstAngles(i,2),'*k','MarkerSize',15,'linewidth',1.5);
    end
    
    xlabel('\theta (degrees)');
    ylabel('\phi (degrees)');
    hold off;
    title(['\Phi^{' angularSpectrumMeth '}(\theta,\phi)   |  markers : sources found ']);
end
end
drawnow;
end

