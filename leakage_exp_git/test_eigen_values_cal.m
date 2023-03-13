clear all;

n = 4;

pp = 0.01;
pt = 0.02;

A1 = pp + pt - (n+1)*2^n*pp * pt;
B1 = 1 - 2*(pp + pt) + 2^(n+2) * pp * pt;
B2 = 2^(n+1) * pp * pt;

fprintf('a1=%f, b1 = %f, b2 = %f\n', A1, B1, B2);
Q = [1-n * A1, 1 - (n-1)*B2 - B1,  1- (n-1)*B2 - B1, 1- (n-1)*B2 - B1, 1- (n-1)*B2 - B1;
    A1,         B1,                 B2,             B2,                 B2;
    A1,         B2,                 B1,             B2,                 B2;
    A1,         B2,                 B2,             B1,                 B2;
    A1,         B2,                 B2,             B2,                 B1         
];

lam_c1 = 1 - 2*(pp + pt) + 2^(n+1) * pp * pt;
lam_c2 = 1 - (n+2) * (pt + pp) +  (n+2)*(n+1)*pp * pt * 2^n;

fprintf('cal_1: %.4f, caL_2:%.4f\n',lam_c1,lam_c2);
[V,D] = eig(Q);
fprintf('%.4f ',diag(D));
