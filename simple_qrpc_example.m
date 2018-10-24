
%SIMPLE_QRPC_EXAMPLE a simple synthetic example demonstrating
% the use of QRPC routines
%
% See also QRPC

% MIT License

% random number generator
iseed = 8675309;
rng(iseed);

% sizes 
m = 100; % number of samples
n = 200; % number of sensor locations
k = 40; % approximate rank
kmax = k; % max number of sensors
kslice = 40; % number of sensors for comparison purposes...
ngamma = 20; % number of relative weightings of fit and cost
gamma_scale = 0.1; % scale of relative weights
eta = 1e-1; % a background noise level

X = randn(m,k)*[randn(k,n-k), eye(k)]+eta*[randn(m,n-k), zeros(m,k)]; 
                                    % data is almost spanned by 
                                    % last k columns
c1 = 1:n; % sensor costs

                % perform sensor placement for different relative weights
gammas = linspace(0,1,ngamma)*gamma_scale;
gammas = gammas(ngamma:-1:1);
errs = zeros(ngamma,1);
errsbig = zeros(ngamma,kmax);
costs = zeros(ngamma,1);
costsbig = zeros(ngamma,kmax);

tic
for i = 1:length(gammas)
    gamma = gammas(i);
    c = c1*gamma;
    if exist('xqrmc_m','file')==3
        [aref,p,tau] = xqrmc_m(X,c,kmax);
        errs(i) = recon_v2(aref,p,X,kslice);
         for j = 1:kmax
            errsbig(i,j) = recon_v2(aref,p,X,j);
        end

    else
        [~,~,p] = qrpc(X,c,kmax);
        errs(i) = recon(X.',X.',p(1:kslice));        
        for j = 1:kmax
            errsbig(i,j) = recon(X.',X.',p(1:j));
        end
    end
    
    costs(i) = sum(c1(p(1:kslice)));     
    costsbig(i,:) = cumsum(c1(p(1:kmax)));

end
toc
%

figure(1)

subplot(2,2,1)
plot(gammas,errs,'x')
title(sprintf('error as gamma increases\n fixed k = %d',kslice))
subplot(2,2,2)
plot(gammas,costs,'x')
title(sprintf('cost as gamma increases\n fixed k = %d',kslice))
subplot(2,2,3)
plot(costs,errs,'x')
title(sprintf('error vs cost\n fixed k = %d',kslice))

subplot(2,2,4)
cmin = min(costsbig(:));
cmax = max(costsbig(:));
ctest = linspace(cmin,cmax);
etest = zeros(size(ctest));
for i = 1:length(ctest)
    etest(i) = min(min(errsbig(costsbig <= ctest(i))));
end

semilogy(ctest,etest,'x')
title(sprintf('error vs cost\n variable k'))