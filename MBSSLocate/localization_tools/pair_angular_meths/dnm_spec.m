function spec = dnm_spec(hatRxx, f, d, tauGrid)

% APR_SPEC Computes the SNR in all directions using ML under a diffuse
% noise model
%
% spec = dnm_spec(hatRxx, f, d, tauGrid)
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
R11 = real(hatRxx(:,:,1,1));
R12 = hatRxx(:,:,1,2);
R21 = hatRxx(:,:,2,1);
R22 = real(hatRxx(:,:,2,2));
c = 343;
SINC = sinc(2*f*d/c);
SINC2 = SINC.^2;

% Initializing the variances
vs = zeros(nbin,nFrames,ngrid);
vb = zeros(nbin,nFrames,ngrid);
for pkInd = 1:ngrid,
    
    % Computing inv(A) = [invA11 invA12; conj(invA11) -invA12]
    EXP = exp(-2*1i*pi*tauGrid(pkInd)*f);
    P = SINC .* EXP;
    invA11 = sqrt(.5)./(1-real(P)).*(1-conj(P));
    invA12 = -(1-P)./(SINC-EXP).*invA11;
    
    % Computing inv(Lambda) = [.5 invL12; 0 invL22]
    DEN = .5./(1-2*real(P)+SINC2);
    invL12 = (SINC2-1).*DEN;
    invL22 = 2*(1-real(P)).*DEN;
    
    % Computing vs and vb without nonnegativity constraint
    ARA1 = repmat(abs(invA11).^2,1,nFrames).*R11 + repmat(abs(invA12).^2,1,nFrames).*R22;
    ARA2 = ARA1 - 2 * real(repmat(invA11.*invA12,1,nFrames).*R21);
    ARA1 = ARA1 + 2 * real(repmat(invA11.*conj(invA12),1,nFrames).*R12);
    vsind = .5*ARA1 + repmat(invL12,1,nFrames).*ARA2;
    vbind = repmat(invL22,1,nFrames).*ARA2;
    
    % Enforcing the nonnegativity constraint (on vs only)
    neg = (vsind < 0) | (vbind < 0);
    vsind(neg) = 0;
    vs(:,:,pkInd) = vsind;
    vb(:,:,pkInd) = vbind;
end
spec = vs./vb;

end