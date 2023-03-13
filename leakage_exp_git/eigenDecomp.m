%output the eigenstates and eigenvalues of the parameterized matrix

clear all;

% syms eps p_l p_s

eps = 0.01;
p_l = 0.02;
p_s = 0.02;

U_1 = [1/2 1/2 0;1/2 1/2 0;0 0 1];

U = kron(U_1, U_1);
E = [
1 0 0 0 0 0 0 0 0;
0 0 0 1 0 0 0 0 0;
0 0 1-eps/2 0 eps/2 0 0 0 0;
0 1 0 0 0 0 0 0 0;
0 0 eps/2 0 1-eps 0 eps/2 0 0;
0 0 0 0 0 1 0 0 0;
0 0 0 0 eps/2 0 1-eps/2 0 0;
0 0 0 0 0 0 0 1 0;
0 0 0 0 0 0 0 0 1
];

Ibb = eye(9);
vec_I1 = [1,1,0,1,1,0,0,0,0];
vec_I2 = [0,0,1,0,0,1,1,1,1];
d1 = 4;
d2 = 5;

Ed = (1 - p_l - p_s)* Ibb + (p_s * vec_I1'/d1 + p_l * vec_I2'/d2)*ones(1,9);

 App_U_1= U * Ed * E;
 
 display(App_U_1);
[states, values] = eig(App_U_1);

values = diag(values);






% Ed = [
% 1,  0,      0,   0,         0,            0,     0,     0,            0;
% 0,  1-p2,  p2,   0,         0,            0,     0,     0,            0;
% 0,  p2,  1-p2,   0,         0,            0,     0,     0,            0;
% 0,  0,      0, 1-p1,        0,            0,    p1,     0,            0;
% 0,  0,      0,  0,      3-p1-p2-p12,     p2,     0,     p1,          p12;
% 0,  0,      0,  0,          p2,     2-p1-p2,     0,     0,           p1;
% 0,  0,      0,  p1,         0,            0,   1-p1,    0,            0;
% 0,  0,      0,  0,          p1,           0,      0, 2-p1-p2,        p2;
% 0,  0,      0,  0,         p12,          p1,      0,    p2,  3-p1-p2-p12
% ];


