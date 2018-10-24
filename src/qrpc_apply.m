function Y = qrpc_apply(QH,X,varargin)
%QRPC_APPLY apply Householder reflectors as returned by QRPC
%
% Input:
%   QH - matrix storing reflectors, represents orthogonal matrix Q
%   X - matrix to multiply
%   varargin{1} - string, determines whether you should
%                 form Q*X or Q'*X. If varargin{1} == 't' then
%                 it's Q'*X. Default is Q*X
% Output:
%   Y - either Q*X or Q'*X
%
% Examples:
%   >> Y = qrpc_apply(QH,X,'t'); % gives Q'*X
%   >> Y = qrpc_apply(QH,X); % gives Q*X
%
% Reference: GW Stewart "Matrix Algorithms, Volume I: Basic 
% Decompositions", Cleve's Corner 
% http://blogs.mathworks.com/cleve/2016/10/03/ ... 
%   householder-reflections-and-the-qr-decomposition/
% (remove space and ellipses)
%
  
% Copyright Travis Askham 2017
% Available freely under the MIT License
%
% See also QRPC

Y = X;
[~,n] = size(QH);

if (nargin > 2)
    T = varargin{1};
else
    T = 'n';
end

if (T == 't' || T == 'T')
    for j = 1:n
        Y = Y - bsxfun(@times,QH(:,j),QH(:,j)'*Y);
    end
else    
    for j = n:-1:1
        Y = Y- bsxfun(@times,QH(:,j),QH(:,j)'*Y);
    end
end

end
