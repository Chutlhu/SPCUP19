function specInst = MBSS_computeAngularSpectrum(functionName,debugMode,aggregationParam,X,f,freqBins)

% Function MBSS_computeAngularSpectrum
%
% This function call the localization method designed by functionName with
% aggregationParam, X and f as static argument + varargin arguments
%
% INPUT (mandatory):
% functionName     : string, the multi-channel angular spectrum method
% debugMode        : 1 x 1, flag to display local angular spectra
% aggregationParam : struct, containning the 2-channel angular spectra
%                    aggregation parameter. These parameters are the following :
%   * pairId       : nPairs x 2, All microphone pair indexes
%   * d            : nPairs x 1, For each pair, distance (in meters)
%                    between microphones
%   * alpha        : nPairs x nGrid : Array of angles for each
%                    microphone pair corresponding to all {azimut, elevation}
%                    to be tested.
%   * alphaSampled : 1 x nPairs cell array, each cell element contains the
%                    uniformly distributed angles to be tested for the
%                    corresponding pair
%   * azimuth      : 1 x nAz, vector of azimuth values
%   * elevation    : 1 x nEl, vector of elevation values
%   * tauGrid      : 1 x nMicPair cell array, each cell element contains the
%                    TDOA corresponding to the alphaSampled for each pair
%   * nGrid        : Size of the global grid
%   * nPairs       : 1 x 1, number of pairs
%   * c            : 1 x 1, speed of sound (m/s)
% X                : nfreq x nfram x N : N multichannel time-frequency transformed signals
%                                                        OR
%                    nfreq x nfram x N x N : spatial covariance matrices in all time-frequency bins
% f                : nfreq x 1, frequency axis in Hertz
% freqBins         : 1 x K containing the index of frequency bins used for the aggregation

% OUTPUT:
% specInst : nGrid x nFrames, angular spectrum for each frame
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

%
switch functionName
    case 'GCC-PHAT'
        specInst = GCC_PHAT(aggregationParam,X,f,freqBins,debugMode);
    case 'GCC-NONLIN'
        specInst = GCC_NONLIN(aggregationParam,X,f,freqBins,debugMode);
    case 'MVDR'
        specInst = MVDR(aggregationParam,X,f,freqBins,debugMode); % X is a hatRxx in this case
    case 'MVDRW'
        specInst = MVDRW(aggregationParam,X,f,freqBins,debugMode); % X is a hatRxx in this case
    case 'DS'
        specInst = DS(aggregationParam,X,f,freqBins,debugMode); % X is a hatRxx in this case
    case 'DSW'
        specInst = DSW(aggregationParam,X,f,freqBins,debugMode); % X is a hatRxx in this case
    case 'DNM'
        specInst = DNM(aggregationParam,X,f,freqBins,debugMode); % X is a hatRxx in this case
end

end

%% MULTICHANNEL ANGULAR SPECTRUM METHODS
function [specInst] = GCC_PHAT(aggregationParam,X,f,freqBins,debugMode)

% Function GCC_PHAT
%
% Compute the GCC-PHAT algorithm for all pairs of microphones.
%
% INPUT:
% aggregationParam: structure (see the description above)
% X :        nfreq x nfram x N , N multichannel time-frequency transformed
%            signals
% f:         1 x nfreq, frequency axis
% freqBins:  1 x K containing the index of frequency bins used for the aggregation
% debugMode: 1 x 1, flag to display local angular spectra
%
% OUTPUT:
% specInst:   nGrid x nFrames, angular spectrum for each frame
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

%% Computing the angular spectrum
[~,nFrames,~] = size(X);
specInst = zeros(aggregationParam.nGrid, nFrames);

for i = 1:aggregationParam.nPairs
    spec = phat_spec(X(freqBins,:,aggregationParam.pairId(i,:)), f(freqBins), aggregationParam.tauGrid{i}); % NV % [freq x fram x local angle for each pair]
    % sum on frequencies
    specSampledgrid = (shiftdim(sum(spec,1)))';
    
    % Order 1 interpolation on the entire grid
    specCurrentPair = interp1q(aggregationParam.alphaSampled{i}', specSampledgrid, aggregationParam.alpha(i,:)');
    
    % Aggregation
    specInst(:,:) = specInst(:,:) + specCurrentPair;
    
    if(debugMode)
        if(i == 1)
            figId = figure;
            colormap(jet); % bleu jaune rouge
        end
        plot_debug(figId,specSampledgrid,specCurrentPair,specInst,aggregationParam.pairId(i,:),aggregationParam.azimuth,aggregationParam.elevation,aggregationParam.tauGrid{i});
    end
end

end

function [specInst] = GCC_NONLIN(aggregationParam,X,f,freqBins,debugMode)

% Function GCC_NONLIN
%
% Compute the GCC-NONLIN algorithm for all pairs of microphones.
%
% INPUT:
% aggregationParam: structure (see the description above)
% X:                nfreq x nfram x N , N multichannel time-frequency
%                   transformed signals
% f:                1 x nfreq, frequency axis
% freqBins:         1 x K containing the index of frequency bins used for
%                   the aggregation
% debugMode:        1 x 1, flag to display local angular spectra
%
% OUTPUT:
% specInst:         1 x nDirection, angular spectrum
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

%% Computing the angular spectrum
% local spectrum
alpha_meth = (10*aggregationParam.c)./(aggregationParam.d*aggregationParam.fs);
[~,nFrames,~] = size(X);
specInst = zeros(aggregationParam.nGrid, nFrames);

for i = 1:aggregationParam.nPairs
    spec = nonlin_spec(X(freqBins,:,aggregationParam.pairId(i,:)), f(freqBins), alpha_meth(i), aggregationParam.tauGrid{i});
    % sum on frequencies
    specSampledgrid = (shiftdim(sum(spec,1)))';
    % Order 1 interpolation on the entire grid
    specCurrentPair = interp1q(aggregationParam.alphaSampled{i}', specSampledgrid, aggregationParam.alpha(i,:)');
    
    %Aggregation
    specInst = specInst + specCurrentPair;
    
    if(debugMode)
        if(i == 1)
            figId = figure;
            colormap(jet); % bleu jaune rouge
        end
        plot_debug(figId,specSampledgrid,specCurrentPair,specInst,aggregationParam.pairId(i,:),aggregationParam.azimuth,aggregationParam.elevation,aggregationParam.tauGrid{i});
    end
end

end

function [specInst] = MVDR(aggregationParam,hatRxx, f,freqBins,debugMode)

% Function MVDR_MULTI
%
% Compute the MVDR_SPEC algorithm for all pairs of microphones.
%
% INPUT:
% aggregationParam: structure (see the description above)
% hatRxx:           nfreq x nfram x N x N , spatial covariance matrices in all
%                   time-frequency bins
% f:                1 x nfreq, frequency axis
% freqBins:         1 x K containing the index of frequency bins used for 
%                   the aggregation
% debugMode:        1 x 1, flag to display local angular spectra
%
% OUTPUT:
% specInst:         1 x nDirection, angular spectrum
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

%% Computing the angular spectrum
[~,nFrames,~,~] = size(hatRxx); % nbin x nFrames x 2 x 2
specInst = zeros(aggregationParam.nGrid, nFrames);

for i = 1:aggregationParam.nPairs
    spec = mvdr_spec(hatRxx(freqBins,:,aggregationParam.pairId(i,:),aggregationParam.pairId(i,:)), f(freqBins), aggregationParam.tauGrid{i}); %
    % sum on frequencies
    specSampledgrid = (shiftdim(sum(spec,1)))';
    
    % Order 1 interpolation on the entire grid
    specCurrentPair = interp1q(aggregationParam.alphaSampled{i}', specSampledgrid, aggregationParam.alpha(i,:)');
    
    %Aggregation
    specInst = specInst + specCurrentPair;
    
    if(debugMode)
        if(i == 1)
            figId = figure;
            colormap(jet); % bleu jaune rouge
        end
        plot_debug(figId,specSampledgrid,specCurrentPair,specInst,aggregationParam.pairId(i,:),aggregationParam.azimuth,aggregationParam.elevation,aggregationParam.tauGrid{i});
    end
end

end

function [specInst] = MVDRW(aggregationParam,hatRxx, f,freqBins,debugMode)

% Function MVDRW_MULTI
%
% Compute the MVDRW_SPEC algorithm for all pairs of microphones.
%
% INPUT:
% aggregationParam: structure (see the description above)
% hatRxx:           nfreq x nfram x N x N , spatial covariance matrices in all
%                   time-frequency bins
% f:                1 x nfreq, frequency axis
% freqBins:         1 x K containing the index of frequency bins used for 
%                   the aggregation
% debugMode:        1 x 1, flag to display local angular spectra
%
% OUTPUT:
% specInst:         1 x nDirection, angular spectrum
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

%% Computing the angular spectrum
% local spectrum
[~,nFrames,~,~] = size(hatRxx); % nbin x nFrames x 2 x 2
specInst = zeros(aggregationParam.nGrid, nFrames);

for i = 1:aggregationParam.nPairs
    spec = mvdrw_spec(hatRxx(freqBins,:,aggregationParam.pairId(i,:),aggregationParam.pairId(i,:)), f(freqBins), aggregationParam.d(i), aggregationParam.tauGrid{i}); %
    % sum on frequencies
    specSampledgrid = (shiftdim(sum(spec,1)))';
    
    % Order 1 interpolation on the entire grid
    specCurrentPair = interp1q(aggregationParam.alphaSampled{i}', specSampledgrid, aggregationParam.alpha(i,:)');
    
    %Aggregation
    specInst = specInst + specCurrentPair;
    if(debugMode)
        if(i == 1)
            figId = figure;
            colormap(jet); % bleu jaune rouge
        end
        plot_debug(figId,specSampledgrid,specCurrentPair,specInst,aggregationParam.pairId(i,:),aggregationParam.azimuth,aggregationParam.elevation,aggregationParam.tauGrid{i});
    end
end

end

function [specInst] = DS(aggregationParam, hatRxx, f,freqBins,debugMode)

% Function DS_MULTI
%
% Compute the DS_SPEC algorithm for all pairs of microphones.
%
% INPUT:
% aggregationParam: structure (see the description above)
% hatRxx:           nfreq x nfram x N x N , spatial covariance matrices in all
%                   time-frequency bins
% f:                1 x nfreq, frequency axis
% freqBins:         1 x K containing the index of frequency bins used for 
%                   the aggregation
% debugMode:        1 x 1, flag to display local angular spectra
%
% OUTPUT:
% specInst:         1 x nDirection, angular spectrum
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

%% Computing the angular spectrum
% local spectrum
[~,nFrames,~,~] = size(hatRxx); % nbin x nFrames x 2 x 2
specInst = zeros(aggregationParam.nGrid, nFrames);

for i = 1:aggregationParam.nPairs
    spec = ds_spec(hatRxx(freqBins,:,aggregationParam.pairId(i,:),aggregationParam.pairId(i,:)), f(freqBins), aggregationParam.tauGrid{i}); %
    % sum on frequencies
    specSampledgrid = (shiftdim(sum(spec,1)))';
    
    % Order 1 interpolation on the entire grid
    specCurrentPair = interp1q(aggregationParam.alphaSampled{i}', specSampledgrid, aggregationParam.alpha(i,:)');
    
    %Aggregation
    specInst = specInst + specCurrentPair;
    
    if(debugMode)
        if(i == 1)
            figId = figure;
            colormap(jet); % bleu jaune rouge
        end
        plot_debug(figId,specSampledgrid,specCurrentPair,specInst,aggregationParam.pairId(i,:),aggregationParam.azimuth,aggregationParam.elevation,aggregationParam.tauGrid{i});
    end
end

end

function [specInst] = DSW(aggregationParam, hatRxx, f,freqBins,debugMode)

% Function DSW_MULTI
%
% Compute the DSW_SPEC algorithm for all pairs of microphones.
%
% INPUT:
% aggregationParam: structure (see the description above)
% hatRxx:           nfreq x nfram x N x N , spatial covariance matrices in all
%                   time-frequency bins
% f:                1 x nfreq, frequency axis
% freqBins:         1 x K containing the index of frequency bins used for 
%                   the aggregation
% debugMode:        1 x 1, flag to display local angular spectra
%
% Outputs :
% specInst:         1 x nDirection, angular spectrum
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

%% Computing the angular spectrum
% local spectrum
[~,nFrames,~,~] = size(hatRxx); % nbin x nFrames x 2 x 2
specInst = zeros(aggregationParam.nGrid, nFrames);

for i = 1:aggregationParam.nPairs
    spec = dsw_spec(hatRxx(freqBins,:,aggregationParam.pairId(i,:),aggregationParam.pairId(i,:)), f(freqBins), aggregationParam.d(i), aggregationParam.tauGrid{i}); %
    % sum on frequencies
    specSampledgrid = (shiftdim(sum(spec,1)))';
    
    % Order 1 interpolation on the entire grid
    specCurrentPair = interp1q(aggregationParam.alphaSampled{i}', specSampledgrid, aggregationParam.alpha(i,:)');
    
    %Aggregation
    specInst = specInst + specCurrentPair;
    
    if(debugMode)
        if(i == 1)
            figId = figure;
            colormap(jet); % bleu jaune rouge
        end
        plot_debug(figId,specSampledgrid,specCurrentPair,specInst,aggregationParam.pairId(i,:),aggregationParam.azimuth,aggregationParam.elevation,aggregationParam.tauGrid{i});
    end
end

end

function [specInst] = DNM(aggregationParam, hatRxx, f,freqBins,debugMode)

% Function DNM_MULTI
%
% Compute the DNM_SPEC algorithm for all pairs of microphones.
%
% INPUT:
% aggregationParam: structure (see the description above)
% hatRxx:           nfreq x nfram x N x N , spatial covariance matrices in all
%                   time-frequency bins
% f:                1 x nfreq, frequency axis
% freqBins:         1 x K containing the index of frequency bins used for 
%                   the aggregation
% debugMode:        1 x 1, flag to display local angular spectra
%
% OUTPUT:
% specInst : 1 x nDirection, angular spectrum
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

%% Computing the angular spectrum
% local spectrum
[~,nFrames,~,~] = size(hatRxx); % nbin x nFrames x 2 x 2
specInst = zeros(aggregationParam.nGrid, nFrames);

for i = 1:aggregationParam.nPairs
    spec = dnm_spec(hatRxx(freqBins,:,aggregationParam.pairId(i,:),aggregationParam.pairId(i,:)) , f(freqBins), aggregationParam.d(i), aggregationParam.tauGrid{i}); %
    % sum on frequencies
    specSampledgrid = (shiftdim(sum(spec,1)))';
    
    % Order 1 interpolation on the entire grid
    specCurrentPair = interp1q(aggregationParam.alphaSampled{i}', specSampledgrid, aggregationParam.alpha(i,:)');
    
    %Aggregation
    specInst = specInst + specCurrentPair;
    
    if(debugMode)
        if(i == 1)
            figId = figure;
            colormap(jet); % bleu jaune rouge
        end
        plot_debug(figId,specSampledgrid,specCurrentPair,specInst,aggregationParam.pairId(i,:),aggregationParam.azimuth,aggregationParam.elevation,aggregationParam.tauGrid{i});
    end
end

end

%% Functions could be used for debug purpose
function plot_debug(figId,specSampledgrid,specCurrentPair,specEntireGrid,pairId,azimuth,elevation,tauGrid)
% This function is used to plot:
% (a) the sampled angular spectrum into the current microphone referential as a function of TDOA
% (b) the angular spectrum interpolated over the entire search grid
% (c) the pooled angular spectrum over pairs

micPairSpecSampledGrid = shiftdim(max(specSampledgrid,[],2));
micPairSpecEntireGrid = shiftdim(max(squeeze(specCurrentPair(:,:)),[],2));
specInstEntireGrid = shiftdim(max(squeeze(specEntireGrid(:,:)),[],2));

iNbThetas = length(azimuth);
iNbPhis = length(elevation);

% (a) mic pair sampled spec
figure(figId);
subplot(1,3,1);
x_axes = tauGrid;
plot(x_axes,micPairSpecSampledGrid);
title({'Angular spectrum in the microphone pair referential' ; ['(microphones ' num2str(pairId(1)) '-' num2str(pairId(2)) ')'] ; 'Sampled'});
xlabel('TDOA (sec)');

% (b) mic pair spec
subplot(1,3,2);
ppfmicPairSpec2D = (reshape(micPairSpecEntireGrid,iNbThetas,iNbPhis))';
imagesc(azimuth,elevation,ppfmicPairSpec2D);
title({'Angular spectrum in the microphone pair referential' ; ['(microphones ' num2str(pairId(1)) '-' num2str(pairId(2)) ')'] ; 'Interpolated over the entire search grid'});
xlabel('Azimuth (°)')
ylabel('Elevataion (°)');

% (c) Aggregated spec
subplot(1,3,3);
ppfInstSpec2D = (reshape(specInstEntireGrid,iNbThetas,iNbPhis))';
imagesc(azimuth,elevation,ppfInstSpec2D);
title('Aggregated angular spectrum over pairs');
xlabel('Azimuth (°)')
ylabel('Elevataion (°)');

drawnow;
end