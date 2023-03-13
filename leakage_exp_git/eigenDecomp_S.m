%%%%here we let pd1 = 0.1 * pd1 in eigenDecomp_General.m
clear all;
%p2--> p2, pd2-->pd2
syms pd1

p1 = 0;
p2 = 0;
pd1 = 0;
pd3 = 0;


p3 = 0;

e_0 = [1 0 0]';
e_1 = [0 1 0]';
e_2 = [0 0 1]';

E= kron(e_0,e_0) * kron(e_0,e_0)' + (1-p2)* kron(e_0,e_1) * kron(e_0,e_1)'...
    +(1-p2-pd1)* kron(e_0,e_2) * kron(e_0,e_2)'...
    +(1-p1)* kron(e_1,e_0) * kron(e_1,e_0)'...
    +(1-p1 - p2-p3-pd1-pd2)* kron(e_1,e_1) * kron(e_1,e_1)'...
    +(1-p1 -p2 -pd3)* kron(e_1,e_2) * kron(e_1,e_2)'...
    +(1-p1 -pd2)* kron(e_2,e_0) * kron(e_2,e_0)'...
    +(1-p1 -p2 - pd3)* kron(e_2,e_1) * kron(e_2,e_1)'...
    +(1-p1 -p2-p3)* kron(e_2,e_2) * kron(e_2,e_2)';

E = E + p1*(kron(e_1,e_0)*kron(e_2,e_0)' + kron(e_2,e_0)*kron(e_1,e_0)'...
   +kron(e_1,e_1)*kron(e_2,e_1)' + kron(e_2,e_1)*kron(e_1,e_1)'...
   +kron(e_1,e_2)*kron(e_2,e_2)' + kron(e_2,e_2)*kron(e_1,e_2)'); %1*<->2*

E = E + p2*(kron(e_0,e_1)*kron(e_0,e_2)' + kron(e_0,e_2)*kron(e_0,e_1)'...
   +kron(e_1,e_1)*kron(e_1,e_2)' + kron(e_1,e_2)*kron(e_1,e_1)'...
   +kron(e_2,e_1)*kron(e_2,e_2)' + kron(e_2,e_2)*kron(e_2,e_1)');%*1<->*2

%11<->22
E = E + p3 * (kron(e_1,e_1)*kron(e_2,e_2)' + kron(e_2,e_2)*kron(e_1,e_1)');

%02<->11
E = E + pd1 * (kron(e_0,e_2)*kron(e_1,e_1)' + kron(e_1,e_1)*kron(e_0,e_2)');

%11<->20
E = E + pd2 * (kron(e_1,e_1)*kron(e_2,e_0)' + kron(e_2,e_0)*kron(e_1,e_1)');

%12<->21
E = E + pd3 * (kron(e_1,e_2)*kron(e_2,e_1)' + kron(e_2,e_1)*kron(e_1,e_2)');

U_1 = [1/2 1/2 0;1/2 1/2 0;0 0 1];
P = kron(U_1, U_1);

UF = P * E;

[states, vals] = eig(UF);


vals = diag(vals);

