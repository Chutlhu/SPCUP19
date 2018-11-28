function Rxx = MBSS_covariance(X)
% Function MBSS_covariance
% This function computes the averaged covariance of X over time axis.
%
% INPUT:
% X: F x N x I, the multichannel time-frequency representation of x.
%  F: number of frequency bins
%  N: number of frames
%  I: number of channels
%
% OUTPUT:
% Rxx: I x I x F, the averaged covariance of X
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

[F,N,I] = size(X);
Rxx = zeros(I,I,F);
for f = 1:F
    for n = 1:N
        Xi = squeeze(X(f,n,:));
        Rxx(:,:,f) = Rxx(:,:,f) +  Xi * Xi';
    end
end
Rxx = Rxx./N;
end