%XORMQR_M MatLab interface for LAPACK QR multiply routine ZUNMQR or DORMQR
%
% Input:
%
% AREF and tau are as returned by xgeqp3_m (should be of same type,
% either real or complex double) These store the necessary information
% for applying the Householder reflectors.
% SIDE = 'L' or 'R', multiply on left or right
% TRANS = 'N' or 'T' or 'C', multiply by Q, Q transpose, or Q conjugate
% transpose
% varargin{1} = k = number of reflectors defining Q (typically k = 
%  min(m,n) where (m,n) are dimensions of AREF)
% B - matrix to be multiplied
%
% Output:
% 
% C as defined below
%
% Example:
%
%   >> C = xormqr_m(SIDE,TRANS,AREF,tau,B,varargin);
%
% C = Q*B or C = B*Q or C = Q'*B or C = B*Q' (determined by SIDE and
% TRANS)
%

% Copyright Travis Askham 2017
% MIT License