function spec = nonlin_spec(X, f, alpha, tauGrid)

% NONLIN_SPEC Computes the nonlinear GCC-PHAT spectrum defined in
% B. Loesch, B. Yang, "Blind source separation based on time-frequency
% sparseness in the presence of spatial aliasing", in 9th Int. Conf. on
% Latent Variable Analysis and Signal Separation (LVA/ICA), pp. 1â€“8, 2010.
%
% spec = nonlin_spec(X, f, alpha, tauGrid)
%
% Inputs:
% X: nbin x nFrames x 2 matrix containing the STFT coefficients of the input
%     signal in all time-frequency bins
% f: nbin x 1 vector containing the center frequency of each frequency bin
%     in Hz
% alpha: nonlinearity parameter
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

spec = zeros(nbin,nFrames,ngrid);
P = X1.*conj(X2);
P = P./abs(P);

temp = ones(1,nFrames);
for pkInd = 1:ngrid,
    EXP = exp(-2*1i*pi*tauGrid(pkInd)*f);
    EXP = EXP(:,temp);
    spec(:,:,pkInd) = 1 - tanh(alpha*sqrt(abs(2-2*real(P.*EXP)))); % RLB : la valeur absolu permet de ne pas se retrouver avec des spec complexe (déjà observé). Je ne connais pas les tenants et aboutissants de ce pb.
end

end