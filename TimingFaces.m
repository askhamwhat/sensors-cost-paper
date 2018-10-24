%TimingFaces this script finds approximate run times for the sensor
% placement routine on the Yale Data set...
% Find sparse sensor locations for reconstruction of eigenfaces,
% incorporating a cost function.

clear all; close all; clc

Z = load('YaleB_32x32.mat');
Z = Z.fea';

[N,L] = size(Z);
n = 32;
r = round(L*.80); % The size of the test set

% Construct the cost function (uncomment your choice)
% %Gaussian
[X1,X2] = meshgrid(-n/2:n/2-1,-n/2:n/2-1);
f2 = exp(-0.05*X1.^2 - 0.05*X2.^2);
f = reshape(f2,N,1);
Gamma = 10^4*[0:0.1:1,1.25,1.5:0.5:4,5:2:19]; % Cost function weightings

kmax = 400;


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

fprintf('running the Gaussian cost function test...\n')
fprintf('dim = %d no. of samples = %d k = %d\n',N,r,kmax)
fprintf('time = %5e no. of trials = %d rate = %5e\n',time,length(Gamma),rate)






