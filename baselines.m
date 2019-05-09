% main function for participants of the SPCUP19 challenge to
% load each recording and to run their algorithm.
%
% Inputs:  N/A (loads data from mat files specified by the
%          variables in PATH)
% Outputs: N/A (saves data to mat files)
%
% Authors: Diego Di Carlo
%          Antoine Deleforge
%
% Notice:  This is a preliminary version as part of the SPCUP19 consultation
%          pre-release. Please report problems and bugs to the author on the
%          PIAZZA platform webpage or on the github page
%

close all
clear
clc


%% DEBUG with EASY DATA
try_with_easy_data = 1;
if try_with_easy_data
    suff = 'e' ;
else
    suff = '';
end

%% METHOD NAME
METHOD_OR_TEAM = [suff 'baseline']; % write here your team/submission name. Please use only one word


%% PATHs
PATH_DATA = [suff 'final/']; % 'challenge
TASK = 'flight';
PATH_AUDIO = [PATH_DATA TASK '_task/audio/'];
PATH_RESULTS = 'results/';

%% EXTERNALS
% add MBSSLocate toolbox to the current matlab session
addpath(genpath('./MBSSLocate/'));

%% GROUND-TRUTH and DEVELOPMENT
% if 'development' is 1, load the ground-truth files and
% perform the evaluation
development = 1;
PATH_GT = ['./../' suff 'final_ground_truth/'];
tollerance = 10;


%% HARD CODED VARIBLEs
%    coord:    x         y         z
micPos = [  0.0420    0.0615   -0.0410;  % mic 1
    -0.0420    0.0615    0.0410;  % mic 2
    -0.0615    0.0420   -0.0410;  % mic 3
    -0.0615   -0.0420    0.0410;  % mic 4
    -0.0420   -0.0615   -0.0410;  % mic 5
    0.0420   -0.0615    0.0410;  % mic 6
    0.0615   -0.0420   -0.0410;  % mic 7
    0.0615    0.0420    0.0410]; % mic 8

%% GROUND TRUTH
if contains(TASK, 'flight')
    offset = 90;
end
if contains(TASK, 'static')
    offset = -90;
end

if development
    % load groud-truth files
    load([PATH_GT TASK '_task/' TASK '_azimuth.mat']);
    load([PATH_GT TASK '_task/' TASK '_elevation.mat']);
    
    if strcmp(TASK,'flight')
        azRef = flight_azimuth;   % azimuth reference
        elRef = flight_elevation; % elevation refernce
    elseif strcmp(TASK,'static')
        azRef = static_azimuth;   % azimuth reference
        elRef = static_elevation; % elevation refernce
    end
    [J, T] = size(azRef); % J audio files x T frames
end


%% BASELINE

% Signal parameters
% -----------------

% Sampling frequency for the algorithm. It can be changed and
% resampling is performed later
fs = 44100;
% if the flight data, than process the signal frame-wise
if strcmp(TASK, 'flight')
    frame_size = 0.500; % 500 [ms] - !hardcode; see Challenge's syllabus
    frame_hop  = 0.250; % 250 [ms] - !hardcode; see Challenge's syllabus
    T = 80;
else
    frame_size = 4; % 4 [s]        - !hardcode; see Challenge's syllabus
    frame_hop  = 0; % no hop size  - !hardcode; see Challenge's syllabus
    T = 1;          % only one frame
end


% MBSS Locate Parameters
% ----------------------

% array parameters
micPos = micPos';     % transpose
isArrayMoving   = 0;  % The microphone array is not moving
subArray        = []; % []: all microphones are used
sceneTimeStamps = []; % Both array and sources are statics => no time stamps
% in the this code, the process is done frame-wise,
% thus MBSSLocate's tracking function is Sdisabled

% localization method
angularSpectrumMeth        = 'GCC-PHAT'; % Local angular spectrum method {'GCC-PHAT' 'GCC-NONLIN' 'MVDR' 'MVDRW' 'DS' 'DSW' 'DNM'}
pooling                    = 'max';      % Pooling method {'max' 'sum'}
applySpecInstNormalization = 0;          % 1: Normalize instantaneous local angular spectra - 0: No normalization
% Search space
azBound                    = [-179 180]; % Azimuth search boundaries ((degree))
elBound                    = [-90   90]; % Elevation search boundaries ((degree))
gridRes                    = 1;          % Resolution (degree) of the global 3D reference system {azimuth,elevation}
alphaRes                   = 5;          % Resolution (degree) of the 2D reference system defined for each microphone pair
% Multiple sources parameters
nsrce                      = 1;          % Number of sources to be detected
minAngle                   = 10;         % Minimum angle between peaks
% Moving sources parameters
blockDuration_sec          = [];         % Block duration in seconds (default []: one block for the whole signal)
blockOverlap_percent       = [];         % Requested block overlap in percent (default []: No overlap) - is internally rounded to suited values
% Wiener filtering
enableWienerFiltering      = 0;          % 1: Process a Wiener filtering step in order to attenuate / emphasize the provided excerpt signal into the mixture signal. 0: Disable Wiener filtering
wienerMode                 = [];         % Wiener filtering mode {'[]' 'Attenuation' 'Emphasis'}
wienerRefSignal            = [];         % Excerpt of the source(s) to be emphasized or attenuated
% Display results
specDisplay                = 0;          % 1: Display angular spectrum found and sources directions found - 0: No display
% Other parameters
speedOfSound               = 343;        % Speed of sound (m.s-1) - typical value: 343 m.s-1 (assuming 20�C in the air at sea level)
fftSize_sec                = [];         % FFT size in seconds (default []: 0.064 sec)
freqRange                  = [];         % Frequency range (pair of values in Hertz) to aggregate the angular spectrum : [] means that all frequencies will be used
% Debug
angularSpectrumDebug       = 0;          % 1: Enable additional plots to debug the angular spectrum aggregation

% Convert algorithm parameters to MBSS structures
% -----------------------------------------------
sMBSSParam = MBSS_InputParam2Struct(angularSpectrumMeth,speedOfSound,fftSize_sec,blockDuration_sec,blockOverlap_percent,pooling,azBound,elBound,gridRes,alphaRes,minAngle,nsrce,fs,applySpecInstNormalization,specDisplay,enableWienerFiltering,wienerMode,freqRange,micPos,isArrayMoving,subArray,sceneTimeStamps,angularSpectrumDebug);

%% FILE-WISE and FRAME-WISE PROCESSING
% variable allocation
J = length(dir([PATH_AUDIO '*.wav']));

azPred = zeros(J, T);
elPred = zeros(J, T);

for j = 1:J
    
    fprintf('Processing audio sample %02i/%02i:\n',j,J)
    
    % Load audio filename
    [wavforms,fs_wav] = audioread([PATH_AUDIO num2str(j) '.wav']);
    % resampling to the given sampling frequency
    [p,q] = rat(fs/fs_wav,0.0001);
    wavforms = resample(wavforms,p,q);
    
    [n_samples, n_chan] = size(wavforms);
    
    % Pick current frame of length frame_size
    for t = 1:T
        
        fprintf('  -- frame: %02i/%02i\n',t,T)
        
        % frame processing
        frame_start = floor(fs*((t-1)*frame_hop))+1;
        frame_end = frame_start + floor(fs*frame_size)-1;
        if frame_end > n_samples
            frame_end = n_samples;
        end
        
        wav_frame = wavforms(frame_start:frame_end,:); % nsampl x nchan
        
        % Run the localization method
        % here you should write your own code
        [azEst, elEst, ~, ~] = ...
            MBSS_locate_spec(wav_frame,wienerRefSignal,sMBSSParam);
        
        % Printing for the development
        if development
            azPred(j,t) = azEst + offset;
            fprintf('Ground-truth:  %1.2f %1.2f\n',azRef(j,t),elRef(j,t));
            elPred(j,t) = elEst;
            fprintf('->Estimation:  %1.2f %1.2f\n\n',azPred(j,t), elPred(j,t));
        end
        
    end
    
    fprintf('\n')
    
end


%% TRAJECTORIES
% plot predicted and ground truth drone trajectories
if development
    
    figure(1)
    if J == 1
        c = ['b'];
    else
        c = parula(J);
    end
    
    for j = 1:J
        if strcmp(TASK, 'flight')
            plot(azRef(j,:),   elRef(j,:), '--', 'color', c(j,:))
            text(azRef(j,8)+3, elRef(j,8)-3, num2str(j))
            hold on
            plot(azPred(j,:), elPred(j,:),'-','color',c(j,:))
            text(azPred(j,8)+3, elPred(j,8)-3, num2str(j))
            title('Ground truth and estimated trajectories')
        else
            scatter( azRef(j), elRef(j),'b','o')
            text( azRef(j)+3, elRef(j)-3, num2str(j))
            hold on
            scatter(azPred(j), elPred(j),'r','x')
            text(azPred(j)+5, elPred(j)-5, num2str(j))
            title('Ground truth and estimated points')
        end
        xlim(azBound)
        ylim(elBound)
        xlabel('Azimuth [degree]')
        ylabel('Elevation [degree]')
    end
end
% fprintf('Press any key to continue...\n')
% pause
% close all

%% EVALUATION
if development == 1
    
    %%% dregon challeng metric
    circ_d = @(pred, ref) acosd(cosd(pred-ref));
    gcirc_d = @(az_pred, az_ref, el_pred, el_ref) ...
        great_circ_dist(1, az_pred, az_ref, el_pred, el_ref);
    
    dregon_score = @(az_pred, az_ref, el_pred, el_ref) ...
        sum(sum(gcirc_d(az_pred, az_ref, el_pred, el_ref) <= tollerance));
    
    dregon_rmse   = @(pred, ref) sqrt(mean(circ_d(pred,ref).^2));
    
    %%% compute metric
    acc   = dregon_score(azPred, azRef, elPred, elRef);
    errAz = dregon_rmse(azPred, azRef);
    errEl = dregon_rmse(elPred, elRef);
    
    fprintf('Dregon Score: %d\n',acc);
    fprintf(' -- RMSE Az.: %1.2f\n',errAz);
    fprintf(' -- RMSE El.: %1.2f\n',errEl);
end

%% PREPARE DATA FOR SUBMISSION
fprintf('Saving data for submission... ')

% check dimenision
if ~all(size(azPred) == [J T] & size(elPred) == [J T])
    error('Dimensions do not match the required ones. They must be (N_files x N_frames)')
end

azimuth = azPred;
elevation = elPred;

mkdir(PATH_RESULTS)
save(['./' PATH_RESULTS 'SPCUP19FINAL_' METHOD_OR_TEAM '_' TASK '.mat'], ...
    'azimuth', 'elevation');

fprintf('done.\n')