% Find sparse sensor locations for reconstruction of cylinder flow,
% incorporating a cost function.

clear all; close all; clc

cyl = load('C:\Users\Emily\Documents\Research\FinalCodes\ALL.mat');

Z = cyl.VORTALL;
m = cyl.m;
n = cyl.n;
N = m*n;
L = size(Z,2);
r = 120; % The size of the training set

% Construct a cost function
f2 = zeros(m,n);
f2(1:99,:) = 1;
f = reshape(f2,N,1);

p = 2:2:20; % The number of sensors
Gamma = [0:0.1:1,1.5:0.5:5,6:1:15]; % The cost function weightings
LCV = 30; % The number of cross validations

% Obtain sensor locations and reconstructions
[Cost,Error,A,stdE,stdC] = CostError(Z,f,Gamma,p,r,LCV);

plot(Error,Cost,'k.')

figure
subplot(1,2,1)
spy(reshape(A(:,1),m,n))
subplot(1,2,2)
spy(reshape(A(:,end),m,n))






