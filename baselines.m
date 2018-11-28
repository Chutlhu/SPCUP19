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

%% PATHs
DATA = 'flight'; % static or flight
PATH_DATA = ['./dev_' DATA '/'];

PATH_AUDIO = [PATH_DATA 'audio/'];
PATH_GT = PATH_DATA;

% add MBSSLocate toolbox to the current matlab session
addpath(genpath('./MBSSLocate/'));

%% FLAGs
% if 'development' is 1, load the ground-truth files and
% perform the evaluation
development = 1;

%% HARD CODED VARIBLEs
%    coord:    x       y       z
micPos = [  0.0615 -0.0420 -0.0410;  % mic 1
            0.0615  0.0420  0.0410;  % mic 2
            0.0420  0.0615 -0.0410;  % mic 3
           -0.0420  0.0615  0.0410;  % mic 4
           -0.0615  0.0420 -0.0410;  % mic 5
           -0.0615 -0.0420  0.0410;  % mic 6
           -0.0420 -0.0615 -0.0410;  % mic 7
            0.0420 -0.0615  0.0410]; % mic 8

%% GROUND TRUTH
if development
    % load groud-truth files
    file = load([PATH_GT 'SPCUP19_dev_' DATA '.mat']);
    
    eval(['gt_azimuth = file.' DATA '_azimuth;'])
    eval(['gt_elevation = file.' DATA '_elevation;'])
    
    [J, T] = size(gt_azimuth); % J audio files x T frames
    
    azRef = gt_azimuth;   % azimuth reference
    elRef = gt_elevation; % elevation refernce
end


%% BASELINE

% Signal parameters
% -----------------

% Sampling frequency for the algorithm. It can be changed and 
% resampling is performed later
fs = 44100; 
% if the flight data, than process the signal frame-wise
if strcmp(DATA, 'flight')
    frame_size = 0.500; % 500 [ms] - !hardcode; see Challenge's syllabus
    frame_hop  = 0.250; % 250 [ms] - !hardcode; see Challenge's syllabus
else
    frame_size = 4; % 4 [s]        - !hardcode; see Challenge's syllabus
    frame_hop  = 0; % no hop size  - !hardcode; see Challenge's syllabus
end


% MBSS Locate Parameters
% ----------------------

% array parameters
alpha = 0;
if strcmp(DATA, 'flight')
    alpha = 180;
end
% rotation matrix to change coordinate references
Rz = [  cosd(alpha) -sind(alpha) 0;
        sind(alpha)  cosd(alpha) 0;
                 0            0  1];
micPos = Rz*micPos';  % apply rotation matrix to coordinates
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
elBound                    = [-90   10]; % Elevation search boundaries ((degree))
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
        [azEst, elEst, block_timestamps,elaps_time] = ...
            MBSS_locate_spec(wav_frame,wienerRefSignal,sMBSSParam);
        
        % Printing for the development
        if development
            azPred(j,t) = azEst;
            fprintf('%1.2f %1.2f\n',azRef(j,t), azPred(j,t));
            elPred(j,t) = elEst;
            fprintf('%1.2f %1.2f\n\n',elRef(j,t), elPred(j,t));
        end
        
    end
    
    fprintf('\n')
    
end


%% TRAJECTORIES
% plot predicted and ground truth drone trajectories
figure(1)
c = parula(J);

if development
    for j = 1:J
        if strcmp(DATA, 'flight')
            plot( azRef(j,:), elRef(j,:),'--','color',c(j,:))
            hold on
            plot(azPred(j,:), elPred(j,:),'-','color',c(j,:))
            text(azPred(j,8)+3, elPred(j,8)-3, num2str(j))
            title('Ground truth and estimated trajectories')
        else
            scatter( azRef(j), elRef(j),'b','o')
            hold on
            scatter(azPred(j), elPred(j),'r','x')
            text(azPred(j)+5, elPred(j)-5, num2str(j))
            title('Ground truth and estimated points')
        end
        xlim([-179 180])
        ylim([-90 0])
        xlabel('Azimuth [degree]')
        ylabel('Elevation [degree]')
    end
end
fprintf('Press any key to continue...\n')
pause
close all
    
%% EVALUATION
if development == 1
    
    %%% dregon challeng metric
    dregon_score  = @(pred, ref) sum(abs(pred-ref)<10,1)/size(ref,1);
    
    %%% compute metric
    azimuth_acc   = dregon_score(azPred, azRef);
    elevation_acc = dregon_score(elPred, elRef);
end

%% PREPARE DATA FOR SUBMISSION
fprintf('Saving data for submission... ')

% check dimenision
if ~all(size(azPred) == [J T] & size(elPred) == [J T])
    error('Dimensions do not match the required ones. They must be (N_files x N_frames)')
end

eval([DATA '_azimuth = azPred;'])
eval([DATA '_elevation = elPred;'])

save(['./SPCUP19_' DATA '.mat'], [DATA '_azimuth'], [DATA '_elevation']);

fprintf('done.\n')