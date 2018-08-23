% Find sparse sensor locations for reconstruction of eigenfaces,
% incorporating a cost function.

clear all; close all; clc

Z = load('C:\Users\Emily\Documents\Research\FinalCodes\YaleB_32x32.mat');
Z = Z.fea';

[N,L] = size(Z);
n = 32; % Images are 32 by 32 pixels
r = round(L*.80); % The size of the test set

% Construct the cost function
f2 = zeros(n,n);

% Gaussian
[X1,X2] = meshgrid(-n/2:n/2-1,-n/2:n/2-1);
f2 = exp(-0.05*X1.^2 - 0.05*X2.^2);
f = reshape(f2,N,1);
Gamma = 10^4*[0:0.1:1,1.25,1.5:0.5:4,5:2:19]; % Cost function weightings

% Block off the center
% f2(11:22,11:22) = 1;
% f = reshape(f2,N,1);
% Gamma = 10^4*[0:0.05:0.5,0.6:0.1:1,1.5:0.5:6];

% Block off left third
% f2(:,1:11) = 1;
% f = reshape(f2,N,1);
% Gamma = 10^4*[0:0.05:0.5,0.6:0.1:1,1.5:0.5:6];


% p = [5,10,15,25,50,75,100:50:400]; % The number of sensors
p = 200;
LCV = 2; % The number of cross validations

% Obtain sensor locations, reconstruction errors, and costs
[Cost,Error,A,stdE,stdC] = CostError(Z,f,Gamma,p,r,LCV);

plot(Cost,Error,'k.')

figure
subplot(1,2,1)
imagesc(reshape(A(:,1),n,n)),colormap hot
subplot(1,2,2)
imagesc(reshape(A(:,end),n,n))















