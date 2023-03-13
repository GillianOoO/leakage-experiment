clear all;

syms p
% 
% p = 0.005;
% %single-qubit case.
U = [1/2 1/2 0;1/2 1/2 0;0 0 1];
E = [1 0 0;0 1-p p; 0 p 1-p];
% 
App_U_1= U * E;
[vectors,eigs] = eig(App_U_1);
eigs = diag(eigs);

% vec2 = orth(vectors)
% lam = 1 - 3 * p/2
% App_U_2= App_U_1 - lam * eye(3)
% [U2,V2] = eig(App_U_2);
% s = det(App_U_2);
% vec = App_U_2 * U(:,3)

% 
% syms p1 p2 p12
% 
% U_1 = [1/2 1/2 0;1/2 1/2 0;0 0 1];
% 
% U = kron(U_1, U_1);
% 
% E = [
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
% 
% App_U_1= U * E;
% % [U, V] = eig(App_U_1)
% 
% lam = 1;
% App_U_2= App_U_1 - lam * eye(9);
% [U2,V2] = eig(App_U_2);
% V2 = diag(V2);
% % display(diag(V2));
% % s = det(App_U_2);
% 
