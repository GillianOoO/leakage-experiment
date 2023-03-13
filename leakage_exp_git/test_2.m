
%syms p q

% mat_p = [1 - 8 * p, 4 * p, 4 * p;
%     4  * p, 1 - 4 * p, 0;
%     4 * p, 0, 1- 4 * p];
% 
% mat_q = [1 - 2 * q, 2 * q, 2 * q;
%     q, 1 - 2 * q, 0;
%     q, 0, 1- 2 * q];
%
% mat = mat_p * mat_q;
n = 2;
pt = 0.01;
pp = 0.02;

A1 = pp + pt - (n+1)*2^n*pp * pt;
B1 = 1 - 2*(pp + pt) + 2^(n+2) * pp * pt;
B2 = 2^(n+1) * pp * pt;
lam = B1-B2;
fprintf('a1=%f, b1 = %f, b2 = %f, lam = %f\n', A1, B1, B2, lam);




mat = [1-n * A1, 1 - (n-1)*B2 - B1,  1- (n-1)*B2 - B1;
    A1,         B1,                 B2;
    A1,         B2,                 B1;    
];

mat_lam = mat - lam * eye(size(mat));

mat_lam * [0;-1;1]

% syms eps1 eps2 eps3
% eps1 = eps1/4;
% eps2 = eps2/4;
% eps3 = 0;
% mat = [1 - eps1 - eps2 , 2 * eps1, 2 * eps2;
%     eps1, 1 - 2 * eps1 - 2 * eps3, 2 * eps3;
%     eps2, 2 * eps3, 1 - 2 * eps2 - 2 * eps3];
% 
[V, D]= eig(mat);
% u = V;
% u(:,1) = V(:,1)/sqrt(dot(V(:,1), V(:,1)));
% for k = 2 : length(u)
%     for j = 1 : k-1
%         u(:,k) = u(:,k) - dot(V(:,k), u(:,j)) * u(:,j);
%     end
%     u(:,k) = u(:,k)/sqrt(dot(u(:,k), u(:,k)));
% end
% 
% display(u);
% D = diag(D);
% for k = 1 : length(D)
%     x = (mat - D(k)* eye(length(D))) * u(:,k);
%     fprintf('%f ', x);
%     fprintf('\n%d\n',k);
% end
