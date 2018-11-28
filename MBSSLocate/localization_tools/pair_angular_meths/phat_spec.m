function spec = phat_spec(X, f, tauGrid)

% PHAT_SPEC Computes the GCC-PHAT spectrum as defined in
% C. Knapp, G. Carter, "The generalized cross-correlation method for
% estimation of time delay", IEEE Transactions on Acoustics, Speech and
% Signal Processing, 24(4):320–327, 1976.
%
% spec = phat_spec(X, f, tauGrid)
%
% Inputs:
% X: nbin x nFrames x 2 matrix containing the STFT coefficients of the input
%     signal in all time-frequency bins
% f: nbin x 1 vector containing the center frequency of each frequency bin
%     in Hz
% tauGrid: 1 x ngrid vector of possible TDOAs in seconds
%
% Output:
% spec: nbin x nFrames x ngrid array of angular spectrum values
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

X1 = X(:,:,1);
X2 = X(:,:,2);

[nbin,nFrames] = size(X1);
ngrid = length(tauGrid);

P = X1.*conj(X2);
P = P./abs(P);

spec = zeros(nbin,nFrames,ngrid);
for pkInd = 1:ngrid
    EXP = repmat(exp(-2*1i*pi*tauGrid(pkInd)*f),1,nFrames);
    %spec(:,:,pkInd) = real(P.*EXP);
    spec(:,:,pkInd) = real(P).*real(EXP) - imag(P).*imag(EXP); % faster !
end

end