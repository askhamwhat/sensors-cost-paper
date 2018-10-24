%TimingSST this routine finds approximate timing info for the SST problem
% Find sparse sensor locations for reconstruction of sea surface temps,
% incorporating a cost function.

clear all; close all; clc

sst = ncread('sst.wkmean.1990-present.nc','sst');
mask = ncread('lsmask.nc','mask');

[m,n] = size(mask);
N = m*n;
L = 1400;
r = 1100; % The size of the training set

for j = 1:L
    sst_mask(:,:,j) = sst(:,:,j).*mask;
    Z(:,j) = reshape(sst_mask(:,:,j),N,1); % The data matrix
end

% Construct a cost function
f2 = ones(m,n);
for j = 3:357
   for jj = 3:177
      if mask(j,jj) == 0
          f2(j-2:j+2,jj-2:jj+2) = 0;
      end
   end
end
f = reshape(f2,N,1);

Zm = zeros(nnz(mask),L);
fm = f(mask(:)>0);
for j = 1:L
    Zm(:,j) = Z(mask(:)>0,j);
end

N = size(Zm,1);

%

Gamma = [0:50:100,125:100:225]; % Cost function weightings
kmax = 300;

X = Zm(:,1:r).'; % take first r samples for timing purposes

tic
for i = 1:length(Gamma)
    c = fm*Gamma(i);
    if exist('xqrmc_m','file')==3
        [AREF,p,TAU] = xqrmc_m(X,c,kmax);
    else
        [QH,R,p] = qrpc(X,c,kmax);
    end
end
time = toc;

%

rate = time/length(Gamma);

fprintf('dim = %d no. of samples = %d k = %d\n',N,r,kmax)
fprintf('time = %5e no. of trials = %d rate = %5e\n',time,length(Gamma),rate)