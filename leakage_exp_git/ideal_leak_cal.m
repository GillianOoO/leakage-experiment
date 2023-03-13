ErrVal = load('ErrVal.txt');

Nq = size(ErrVal,2);
K_lower = 3.864e-6;
K_upper = 0.0005548;
K = 0.0002793;

lambda_lower = exp(-K_upper);
lambda_upper = exp(-K_lower);
lambda = exp(-K);

L = sum(ErrVal(1,:))/2^Nq;

S = sum(ErrVal(2,:))/(3^Nq - 2^Nq);

p = (1-lambda)/3;

Lest = 2 * p;
Sest = 4*p* 8/(81-16);