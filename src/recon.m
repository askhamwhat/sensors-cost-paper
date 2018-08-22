function Error = recon(Psi,Y,ind,N)

theta = Psi(ind,:);
theta_inv = pinv(theta);

for jj = 1:size(Y,2)
    ytrue = Y(:,jj);
    ysparse = ytrue(ind);
    
    % Reconstruct the image
    s = theta_inv*ysparse;
    yrecon = Psi*s;
    
    % Calculate the error
%     E(jj) = norm(abs(yrecon - ytrue))/sqrt(N);
    E(jj) = norm(abs(yrecon - ytrue))/norm(ytrue);
end

Error = mean(E);
