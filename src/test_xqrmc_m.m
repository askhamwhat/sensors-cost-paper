

m = 500;
n= 6000;
A = randn(m,n);
B = A;
c= rand(n,1)*10;

tic; [AREF,JPVT,TAU] = xqrmc_m(A,c); toc
tic; [Q,R1,p1] = qr(A,'vector'); toc

tic; [QH,R,p] = qrpc(A,c); toc

norm(double(JPVT).' - p)
