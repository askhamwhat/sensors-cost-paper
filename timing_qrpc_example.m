
%TIMING_QRPC_EXAMPLE a synthetic example demonstrating the computational 
% complexity of QRPC, or, if available, XQRMC_M (the fast mex version)
%
% This takes a while to run!
%
% See also QRPC XQRMC_M

% MIT License

% random number generator
iseed = 8675309;
rng(iseed);

% parameters
ngamma = 2; % number of relative weightings of fit and cost
gamma_scale = 0.1; % scale of relative weights
eta = 1e-1; % a background noise level

% sizes 

ms = 1000:1000:5000;
ns = 1000:1000:5000;
kpercent = 4:4:20;

times = zeros(length(ms),length(ns),length(kpercent));

for ii = 1:length(ms)
    ii
    for jj = 1:length(ns)
        for kk = 1:length(kpercent)
            m = ms(ii); % number of samples
            n = ns(jj); % number of sensor locations
            k = (m*kpercent(kk))/100; % approximate rank
            X = randn(m,k)*[randn(k,n-k), eye(k)]+eta*[randn(m,n-k), zeros(m,k)]; 
                                                % data is almost spanned by 
                                                % last k columns
            c1 = 1:n; % sensor costs

            % perform sensor placement for different relative weights
            gammas = linspace(0,1,ngamma)*gamma_scale;
            gammas = gammas(ngamma:-1:1);

            tic;
            for i = 1:length(gammas)
                gamma = gammas(i);
                c = c1*gamma;
                if exist('xqrmc_m','file')==3
                    [~,p,~] = xqrmc_m(X,c,k);
                    %Ind(1:p)
                else
                    [~,~,p] = qrpc(X,c,k);
                    %Ind(1:p)
                end
            end
            times(ii,jj,kk) = toc;
        end
    end
end
%

figure(1)

mnk = zeros(length(ms),length(ns),length(kpercent));

for i = 1:length(ms)
    for j = 1:length(ns)
        for k = 1:length(kpercent)
            mnk(i,j,k) = ms(i)*ns(j)*(kpercent(k)*ms(i))/100;
        end
    end
end

hold off
loglog(mnk(:),times(:)/ngamma,'bo')
hold on
loglog(mnk(:),mnk(:)*times(end,end)/ngamma/mnk(end,end),'--r')