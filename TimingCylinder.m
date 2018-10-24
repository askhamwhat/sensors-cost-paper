%TimingCylinder this file gets a typical time for running 
% qrpc or, if available, the faster xqrmc_m (mex file based)
% on the cylinder data.
%
% Should be run twice to eliminate JIT compiling effects

% Find sparse sensor locations for reconstruction of cylinder flow,
% incorporating a cost function.

clear all; close all; clc

cyl = load('CYLINDER_ALL.mat');

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

kmax = 20; % maximum number of sensors
Gamma = [0:0.1:1,1.5:0.5:5,6:1:15]; % The cost function weightings

X = Z(:,1:r).'; % take first r samples for timing purposes

tic
for i = 1:length(Gamma)
    c = f*Gamma(i);
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
