function X = MBSS_wiener(X,Rnn,wienerMode)
% Function MBSS_wiener
% This function filters X in order to attenuate or emphaze a source
% signal represented by its covariance matrix Rnn.
%
% INPUT:
%   X: F x N x I, the multichannel time-frequency representation of mixture x. 
%       F: number of frequency bins
%       N: number of frames
%       I: number of channels
% Rnn: I x I x F, the covariance matrix of the signal excerpt to attenuate
% or emphaze. 
% wienerMode: string, wiener enhencement mode for the localization.
% This apply a "pseudo" wiener filter depending of the chosen mode. 
% Available modes: 'Attenuation' or 'Emphasis'
% 
% OUTPUT:
% X: I x I x F, the averaged covariance of X
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

% Compute the signal covariance matrix
Rxx = MBSS_covariance(X);
[F,N,~] = size(X);

switch(wienerMode)
    case 'Attenuation'
        for f = 1:F
            sigma_f = (Rxx(:,:,f) - Rnn(:,:,f))*inv(Rnn(:,:,f)); 
            for n = 1:N
                X(f,n,:) = sigma_f*squeeze(X(f,n,:));
            end
        end
        
    case 'Emphasis'     
        for f = 1:F
            [U,S,V] = svd(Rxx(:,:,f) - Rnn(:,:,f));
            S(S<=0) = 0;
            S = diag(1./diag(S));
            Rss_inv = V*S*U';
            sigma_f = Rnn(:,:,f)*Rss_inv;
            for n = 1:N
                X(f,n,:) = sigma_f*squeeze(X(f,n,:));
            end
        end
        
    otherwise
        error('Unknown method');
end

end