function [Cx,sampleStart,sampleEnd]=MBSS_qstft_multi(x,fs,wlen,lf,lt)

% File MBSS_qstft_multi.m
% Quadratic linear-scale time-frequency transform based on the local
% covariance of a STFT with sine windows
%
% [Cx,f]=MBSS_qstft_multi(x,fs,wlen,lf,lt)
%
% Inputs:
% x: nsampl x nchan vector containing a multichannel signal
% fs: sampling frequency in Hz
% wlen: length of the STFT window (must be a power of 2)
% lf: half-width of the frequency neighborhood for the computation of
% empirical covariance
% lt: half-width of the time neighborhood for the computation of empirical
% covariance
%
% Output:
% Cx: nchan x nchan x nbin x nfram matrix containing the spatial covariance
% matrices of the input signal in all time-frequency bins
% startSample: nfram x 1 , start sample of each frame
% endSample: nfram x 1, last sample of each frame
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2008-2011 Emmanuel Vincent
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Errors and warnings %%%
if nargin<3, error('Not enough input arguments.'); end
[nsampl,nchan]=size(x);
if nchan>nsampl, error('The input signal must be in columns.'); end
if nargin<4, lf=2; end
if nargin<5, lt=2; end

%%% STFT %%%
[X,sampleStart,sampleEnd]=MBSS_stft_multi(x.',wlen);
[nbin,nfram,nchan]=size(X);

%%% Computation of local covariances for each pair of microphone %%%
winf=hanning(2*lf-1);
wint=hanning(2*lt-1).';

% winf=hanning(2*8-1);
% wint=hanning(2*2-1).';
Cx = zeros(nchan,nchan,nbin,nfram);
%Cx = complex(double(zeros(nchan,nchan,nbin,nfram)));

pairId = nchoosek(1:nchan,2);
[nPairs,~] = size(pairId);

for f=1:nbin,
    for t=1:nfram,
        indf=max(1,f-lf+1):min(nbin,f+lf-1);
        indt=max(1,t-lt+1):min(nfram,t+lt-1);
        nind=length(indf)*length(indt);
        wei=ones(nchan,1)*reshape(winf(indf-f+lf)*wint(indt-t+lt),1,nind);
        XX=reshape(X(indf,indt,:),nind,nchan).';
        local_Cx = (XX.*wei)*XX'/sum(wei(1,:));
        for idPair = 1:nPairs
            Cx(pairId(idPair,:),pairId(idPair,:),f,t) = local_Cx(pairId(idPair,:),pairId(idPair,:));
        end
    end
end

f=fs/wlen*(0:nbin-1).';

return;