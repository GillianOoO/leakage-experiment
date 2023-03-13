clear all;

%syms p1 p2 p3 p4 p5

syms 
p1 = 0;
p2 = 0;
p3 = 0;
p4 = 0;
p5 = 0;

e_0 = [1 0 0]';
e_1 = [0 1 0]';
e_2 = [0 0 1]';

E= kron(e_0,e_0) * kron(e_0,e_0)' + (1-p2)* kron(e_0,e_1) * kron(e_0,e_1)'...
    +(1-p2-p4)* kron(e_0,e_2) * kron(e_0,e_2)'...
    +(1-p1)* kron(e_1,e_0) * kron(e_1,e_0)'...
    +(1-p1 - p2-p3-p4-p5)* kron(e_1,e_1) * kron(e_1,e_1)'...
    +(1-p1 -p2)* kron(e_1,e_2) * kron(e_1,e_2)'...
    +(1-p1 -p5)* kron(e_2,e_0) * kron(e_2,e_0)'...
    +(1-p1 -p2)* kron(e_2,e_1) * kron(e_2,e_1)'...
    +(1-p1 -p2-p3)* kron(e_2,e_2) * kron(e_2,e_2)';

E = E + p1*(kron(e_1,e_0)*kron(e_2,e_0)' + kron(e_2,e_0)*kron(e_1,e_0)'...
   +kron(e_1,e_1)*kron(e_2,e_1)' + kron(e_2,e_1)*kron(e_1,e_1)'...
   +kron(e_1,e_2)*kron(e_2,e_2)' + kron(e_2,e_2)*kron(e_1,e_2)');

E = E + p2*(kron(e_0,e_1)*kron(e_0,e_2)' + kron(e_0,e_2)*kron(e_0,e_1)'...
   +kron(e_1,e_1)*kron(e_1,e_2)' + kron(e_1,e_2)*kron(e_1,e_1)'...
   +kron(e_2,e_1)*kron(e_2,e_2)' + kron(e_2,e_2)*kron(e_2,e_1)');

E = E + p3 * (kron(e_1,e_1)*kron(e_2,e_2)' + kron(e_2,e_2)*kron(e_1,e_1)');
E = E + p4 * (kron(e_0,e_2)*kron(e_1,e_1)' + kron(e_1,e_1)*kron(e_0,e_2)');
E = E + p5 * (kron(e_1,e_1)*kron(e_2,e_0)' + kron(e_2,e_0)*kron(e_1,e_1)');


U_1 = [1/2 1/2 0;1/2 1/2 0;0 0 1];
P = kron(U_1, U_1);

UF = P * E;

[states, vals] = eig(UF);

vals = diag(vals);


