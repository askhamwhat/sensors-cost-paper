function [Cost,Error,A,stdC,stdE] = CostError(Z,f,Gamma2,p2,r,LCV)
% Calculate sensor placement cost and error with a cost function.
% 
% Input:
% 
%   Z - the data set, consisting of L column vectors of length N
%   f - the cost function, a vector of length N
%   Gamma2 - cost function weightings
%   p2 - the number of sensors
%   r - the size of the training set
%   LCV - the number of cross validations
% 
% Output:
% 
%   Cost - the matrix of costs. Cost function weighting varies row by row,
%       number of sensors varies across columns.
%   Error - the matrix of errors, with the same structure as Cost
%   A - the average sensor placements. Consists of vectors of length N,
%       cost function weighting varies across the second dimension, number
%       of sensors varies across the third.
%   stdC - the standard deviation of the cost over the cross validations,
%       with the same structure as Cost
%   stdE - the standard deviation of the error over the cross validations,
%       with the same structure as Cost
% 
% Example:
%
%   >> [Cost,Error,A,stdC,stdE] = CostError(Z,f,Gamma2,p2,r,LCV);
%
% Copyright Emily Clark 2018
% Available freely under the MIT License

[N,L] = size(Z);
Lgamma = length(Gamma2);
Lp = length(p2);

Error = zeros(Lgamma,Lp); Cost = zeros(Lgamma,Lp);
stdE = zeros(Lgamma,Lp); stdC = zeros(Lgamma,Lp);
A = zeros(N,Lgamma,Lp);

for j = 1:Lp % Loop over number of sensors
    tic
    p = p2(j);
    
    for k = 1:Lgamma % Loop over cost function weighting
        Gamma = Gamma2(k);
        
        Error2 = zeros(1,LCV); Cost2 = zeros(1,LCV); A2 = zeros(N,LCV);
        for iter = 1:LCV % Loop for cross validation
            q = randperm(L);
            X = Z(:,q(1:r)); % The interpolative training set
            Y = Z(:,q(r+1:end)); % The interpolative test set
            
            Arand = randn(r,2*p); % Randomized rank reduction
            Psi = X*Arand;
            
            % Obtain sensor locations
            [~,~,Ind] = qrpc(Psi',Gamma*f);
            Ind = Ind(1:p);
            
            % Obtain cost, reconstruction error, and sensor placements
            Error2(iter) = recon(Psi,Y,Ind);
            Cost2(iter) = sum(f(Ind));
            A2(Ind,iter) = 1;
        end
        
        Error(k,j) = mean(Error2);
        Cost(k,j) = mean(Cost2);
        
        stdE(k,j) = std(Error2);
        stdC(k,j) = std(Cost2);

        A(:,k,j) = mean(A2,2);
        
    end
    toc
end


