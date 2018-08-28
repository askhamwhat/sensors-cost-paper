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

Error = norm(Y - (Psi/Psi(ind,:))*Y(ind,:),'fro')/norm(Y,'fro');
