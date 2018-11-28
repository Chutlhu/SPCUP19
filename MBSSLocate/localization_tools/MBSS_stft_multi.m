function [X,startSample,endSample]=MBSS_stft_multi(x,wlen)

% File MBSS_stft_multi.m
% Multichannel short-time Fourier transform (STFT) using
% half-overlapping sine windows.
%
% X=MBSS_stft_multi(x,wlen)
%
% Inputs:
% x: nchan x nsampl matrix containing nchan time-domain mixture signals
% with nsampl samples
% wlen: window length (default: 1024 samples or 64ms at 16 kHz, which is
% optimal for speech source separation via binary time-frequency masking)
%
% Output:
% X: nbin x nfram x nchan matrix containing the STFT coefficients with nbin
% frequency bins and nfram time frames
% startSample: nfram x 1 , start sample of each frame
% endSample: nfram x 1, last sample of each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2008 Emmanuel Vincent
% Copyright 2018 Ewen Camberlein and Romain Lebarbenchon
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Difference with 2008 copyrigth:
% - start and last sample for each frame is provided as output
% - signal is truncated instead of using zero padding to respect window
%   length used (wlen)



%%% Errors and warnings %%%
if nargin<1, error('Not enough input arguments.'); end
if nargin<2, wlen=1024; end
[nchan,nsampl]=size(x);
if nchan>nsampl, error('The signals must be within rows.'); end
if wlen~=4*floor(wlen/4), error('The window length must be a multiple of 4.'); end
if mod(nsampl,wlen/2) ~= 0,  nsampl = floor(nsampl/wlen*2)*wlen/2; x = x(:,1:nsampl); fprintf('Input signal truncated at %d samples with respect to window length used\n',nsampl); end
    
%%% Computing STFT coefficients %%%
% Defining sine window
win=sin((.5:wlen-.5)/wlen*pi).';
nfram=nsampl/wlen*2-1;
nbin=wlen/2+1;

% Start and end sample id of each frames
startSample = ((0:nfram-1).*wlen/2 + 1).';
endSample = ((0:nfram-1).*wlen/2+wlen).';

X = complex(double(zeros(nbin,nfram,nchan)));
for i=1:nchan,
    for t=0:nfram-1,
        % Framing
        frame=x(i,t*wlen/2+1:t*wlen/2+wlen).'.*win;
        % FFT
        fframe=fft(frame);
        X(:,t+1,i)=fframe(1:nbin);
    end
end

return;