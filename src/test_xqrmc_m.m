%TEST_XQRMC_M test whether the XQRMC_M file compiled correctly...

m = 500;
n= 6000;
A = randn(m,n);
c= rand(n,1)*10;

k = 200;

tic; [AREF,JPVT,TAU] = xqrmc_m(A,c,k); toc
tic; [Q,R1,p1] = qr(A,'vector'); toc

tic; [QH,R,p] = qrpc(A,c,k); toc

norm(double(JPVT(1:k)).' - p(1:k))

if (exist('xormqr_m','file')==3)
    qhat = xormqr_m('L','N',AREF(:,1:k),TAU(1:k),eye(m,k));
    rhat = triu(AREF(1:k,:));

    %

    norm(A(:,JPVT(1:k))-qhat*rhat(:,1:k),'fro')/norm(A,'fro')
end