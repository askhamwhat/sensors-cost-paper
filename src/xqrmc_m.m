%XQRMC_M (partial) pivoted QR with a cost function
%
% This routine returns a LAPACK-style description of 
% a partial QR decomposition of the input matrix. The pivots 
% are chosen to maximize the difference between the vector norm 
% and cost (given in a vector input) of the remaining columns 
% at each step, i.e. the greedy maximizer of the log |det| minus the 
% cost.
%
% Input: 
%   X - matrix to decompose
%   c - vector of costs associated to each column (length(c) == size(X,2))
%   k - the size of the partial decomposition
%
% Output:
%   AREF, p, TAU - LAPACK style description of a partial QR 
%                  decomposition
%
%  Q*R ~ X(:,p(1:k)), where R = triu(AREF(:,1:k)) and 
%  Q is stored in the lower triangle of AREF with some 
%  info in TAU. 
%  Q can be recovered using the xormqr_m routine
%
% Example: 
%   >> [AREF,p,TAU] = xqrmc_m(X,c,k);
%   >> Qhat = xormqr_m('L','N',AREF(:,1:k),TAU(1:k),eye(m,k));
%   >> Rhat = triu(AREF(1:k,:));
%   >> norm(X(:,p(1:k))-qhat*rhat(:,1:k),'fro') % should be small
%
% see also XORMQR_M


% MEX File function