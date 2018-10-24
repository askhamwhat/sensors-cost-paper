function [QH,R,p] = qrpc(A,c,varargin)
%QRPC QR decomposition with pivoting and a cost function
%
% Input:
%
%   A - matrix to decompose
%   c - vector of costs (one cost per column of A)
%
% Output:
%
%   QH - stores the Householder reflectors which define 
%        Q in the QR decomposition. To apply Q to a vector
%        see QRPC_APPLY
%   R - upper triangular matrix
%   p - vector of pivots
%
% As noted, the entries of QH define an orthogonal matrix Q (which
% can be accessed via QRPC_APPLY). In the end, we have that 
%
% A(:,p) = Q*R
%
% Reference: GW Stewart "Matrix Algorithms, Volume I: Basic 
% Decompositions", Cleve's Corner 
% http://blogs.mathworks.com/cleve/2016/10/03/ ... 
%   householder-reflections-and-the-qr-decomposition/
% (remove space and ellipses)
%
% Examples:
%
%   >> [QH,R,p] = qrpc(A,c);
%
% See also QRPC_APPLY

% Copyright Travis Askham 2017
% Available freely under the MIT License


[m,n] = size(A);

if (length(c) ~= n)
    error('vector of costs has wrong dimensions')
end

% initialize

QH = zeros(m,n);
R = A;
p = 1:n;

k = min(m,n);
if nargin > 2
    k = varargin{1};
end

for j = 1:k
    [u,ipiv] = qrpc_reflector(R(j:m,j:n),c(p(j:n)));
    % track column pivots
    ipiv = ipiv+j-1;
    itemp = p(j);
    p(j) = p(ipiv);
    p(ipiv) = itemp;
    % switch columns
    temp = R(:,j);
    R(:,j) = R(:,ipiv);
    R(:,ipiv) = temp;
    % apply reflector
    QH(j:m,j) = u;
    R(j:m,j:n) = R(j:m,j:n) - bsxfun(@times,u,u'*R(j:m,j:n));
    R(j+1:m,j) = 0;
end

end

function [u,ipiv] = qrpc_reflector(r,c)
%QRPC_REFLECTOR best reflector with column pivoting and a cost function
% 
% This function generates a Householder reflector
% The pivoting is biased by a cost function, i.e.
% the pivot is chosen as the argmax of norm(r(:,i))-c(i)

% size of each column
dlens = sqrt(sum(r.^2,1));

% choose pivot
[~,ipiv] = max(dlens(:)-c(:));

dlen = dlens(ipiv);

if (dlen > 0.0)
    u = r(:,ipiv)/dlen;
    u(1) = u(1) + sign(u(1)) + (u(1) == 0);
    u = u/sqrt(abs(u(1)));
else
    u = r(:,ipiv);
    u(1) = sqrt(2.0);
end

end