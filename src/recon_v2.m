function Error = recon_v2(AREF,p,Y,k)
%RECON_V2 an efficient method for computing reconstruction
% this routine makes use of the partial QR information in 
% AREF and p. Should follow a call to xqrmc_m
%
% Input:
%
%   AREF, p - as returned by xqrmc_m called on an m x n matrix X
%                  whose rows are samples of a system
%   Y - l x n matrix, test set of samples to reconstruct based on first
%        k entries of p
%   k - number of sensor locations
%
% Output:
%
%   Error - the relative reconstruction error for the test set
%
% Example:
%
%   >> [AREF,p,TAU] = [AREF,p,TAU] = xqrmc_m(X,c,kmax);
%   >> Error = recon_v2(AREF,p,Y,k);
%

%
% Copyright Emily Clark and Travis Askham 2017
% Available freely under the MIT License

rel = norm(Y,'fro');

R1 = triu(AREF(1:k,1:k));
R2 = AREF(1:k,k+1:end);

opts.UT = true;
T1 = linsolve(R1,R2,opts);

[m,n] = size(AREF);

T = zeros(k,n);
T(:,p(k+1:end)) = T1;
T(:,p(1:k)) = eye(k);

Error = norm(Y - Y(:,p(1:k))*T,'fro')/rel;
