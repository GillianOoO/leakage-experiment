%Ecal(rho) = p_0 rho + p1 x1 rho x1 + p2 x2 rho x2 + pd1 xd1 rho xd1 + pd2
%xd2 rho xd2 + pd3 xd3 rho xd3
clear all;

syms pa2 pa4

pa1 = 0;
pa3 = 0;
% pa4 = 0.000001;
% pa2 = 0.000001;
pd2 = 0;
pd1 = 0;
p1 = 0;
p2 = 0;
p3 = 0;
%   pd1 = 0.001; 
% pd2 = 0.0005;
%   pd3 = 0.0001;
pd3 = 0;


% 2*pd3 = 0;
% pd3 = 0;

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

Ea1 = [1 - pa1, 0, pa1;
    0, 1 - pa2, pa2;
    pa1, pa2, 1 - pa1 - pa2];

Ea2 = [1 - pa3, 0, pa3;
    0, 1 - pa4, pa4;
    pa3, pa4, 1 - pa3 - pa4];


UF = P * kron(Ea1, Ea2) * E;

out_1 = 1 - 3/8*(pd1 + pd2) - pd3/2;
out_2 = sqrt(7 * (pd1 - pd2)^2 + 2 * (pd1 - 2 * pd3)^2 + 2 * (pd2 - 2 * pd3)^2)/8;

%fprintf('%.7f %.7f\n', out_1, out_2);

[states, vals] = eig(UF);


vals = diag(vals);
% fprintf('%.7f\t',vals);
% indices = [1,2,3,6,7,8,9];
% state_1 = states(:,indices);
