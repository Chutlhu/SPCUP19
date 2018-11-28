function spec = mvdrw_spec(hatRxx, f, d, tauGrid)

% MVDRW_SPEC Computes the SNR in all directions using the MVDR beamformer
% and frequency weighting
%
% spec = mvdrw_spec(hatRxx, f, d, tauGrid)
%
% Inputs:
% hatRxx : nbin x nFrames x 2 x 2 array containing the spatial covariance
%     matrices of the input signal in all time-frequency bins
% f: nbin x 1 vector containing the center frequency of each frequency bin
%     in Hz
% d: microphone spacing in meters
% tauGrid: 1 x ngrid vector of possible TDOAs in seconds
%
% Output:
% spec: nbin x nFrames x ngrid array of SNR values
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2010-2011 Charles Blandin and Emmanuel Vincent
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)
% If you find it useful, please cite the following reference:
% Charles Blandin, Emmanuel Vincent and Alexey Ozerov, "Multi-source TDOA
% estimation in reverberant audio using angular spectra and clustering",
% Signal Processing 92, pp. 1950-1960, 2012.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nbin,nFrames] = size(hatRxx(:,:,1,1));
ngrid = length(tauGrid);
R11 = hatRxx(:,:,1,1);
R12 = hatRxx(:,:,1,2);
R21 = hatRxx(:,:,2,1);
R22 = hatRxx(:,:,2,2);
TR = real(R11 + R22);
c = 343;
SINC = sinc(2*f*d/c);

SNR = zeros(nbin,nFrames,ngrid);
for pkInd=1:length(tauGrid),
    EXP = repmat(exp(-2*1i*pi*tauGrid(pkInd)*f),1,nFrames);
    NUM = real(R11.*R22 - R12.*R21)./(TR - 2*real(R12.*EXP));
    SNR(:,:,pkInd) = repmat(-(1+SINC)/2,1,nFrames) + repmat((1-SINC)/2,1,nFrames).*NUM./(.5*TR-NUM);
end
spec = SNR;

end