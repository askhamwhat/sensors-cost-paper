function Error = recon(Psi,Y,ind)
% Calculate the reconstruction error for the test set Y given sensors ind
% and basis functions Psi.
%
% Input:
%
%   Psi - the n by r matrix of basis functions
%   Y - the n by m matrix of test set snapshots
%   ind - the 1 by p vector of sensor indices
%
% Output:
%
%   Error - the average reconstruction error for the test set
%
% Example:
%
%   >> Error = recon(Psi,Y,ind);
%
% Copyright Emily Clark 2017
% Available freely under the MIT License

[n,m] = size(Y);

theta = Psi(ind,:);
theta_inv = pinv(theta);

E = zeros(1,m);
for j = 1:m
    ytrue = Y(:,j);
    ysparse = ytrue(ind);
    
    % Reconstruct the snapshot
    s = theta_inv*ysparse;
    yrecon = Psi*s;
    
    % Calculate the error
    E(j) = norm(abs(yrecon - ytrue))/norm(ytrue);
end

Error = mean(E);
