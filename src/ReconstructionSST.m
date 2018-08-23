% Find sparse sensor locations for reconstruction of sea surface temps,
% incorporating a cost function.

clear all; close all; clc

sst = ncread(['C:\Users\Emily\Documents\Research\FinalCodes\'...
    'sst.wkmean.1990-present.nc'],'sst');
mask = ncread(['C:\Users\Emily\Documents\Research\FinalCodes\'...
    'lsmask.nc'],'mask');

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

Gamma = [0:10:100,125:25:225]; % Cost function weightings
p = [5,10,25,50,75,100:50:300]; % The number of sensors
LCV = 10; % The number of cross validations

% Obtain sensor locations and reconstructions
[Cost,Error,A,stdE,stdC] = CostError(Z,f,Gamma,p,r,LCV);

plot(Error,Cost,'k.')

figure
subplot(1,2,1)
spy(reshape(A(:,1),m,n)')
subplot(1,2,2)
spy(reshape(A(:,end),m,n)')









