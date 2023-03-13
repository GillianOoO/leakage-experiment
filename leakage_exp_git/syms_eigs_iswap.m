function mat = syms_eigs_iswap()

P = [
1/2 1/2 0;
1/2 1/2 0;
0 0 1
];

P = kron(P, P);
syms err
% err = 0.001;

E = [
1, 0, 0,        0, 0,       0, 0,       0, 0;
0, 0, 0,        1, 0,       0, 0,       0, 0;
0, 0, 1-err/2,  0, err/2,   0, 0,       0, 0;
0, 1, 0,        0, 0,       0, 0,       0, 0;
0, 0, err/2,    0, 1-err,   0, err/2,   0, 0;
0, 0, 0,        0, 0,       1, 0,       0, 0;
0, 0, 0,        0, err/2,   0, 1-err/2, 0, 0;
0, 0, 0,        0, 0,       0, 0,       1, 0;
0, 0, 0,        0, 0,       0, 0,       0, 1
];

% E_two = [
% 1, 0, 0,        0, 0,       0, 0,       0, 0;
% 0, 1, 0,        0, 0,       0, 0,       0, 0;
% 0, 0, 1-err/2,  0, err/2,   0, 0,       0, 0;
% 0, 0, 0,        1, 0,       0, 0,       0, 0;
% 0, 0, err/2,    0, 1-err,   0, err/2,   0, 0;
% 0, 0, 0,        0, 0,       1, 0,       0, 0;
% 0, 0, 0,        0, err/2,   0, 1-err/2, 0, 0;
% 0, 0, 0,        0, 0,       0, 0,       1, 0;
% 0, 0, 0,        0, 0,       0, 0,       0, 1
% ];

mat = P * E;

% file = 'mat.txt';
% fp = fopen(file, 'w+');
% for i = 1 : size(mat, 2)
%    for j = 1 : size(mat, 1)
%       fprintf(fp, '%', mat(i,j)); 
%    end
% end

save('mat.mat', 'mat');
% save('mat.txt', 'mat', 'sym');
eigs = simplify(eig(mat));


lambda(1) = eigs(9);
lambda(2) = eigs(8);
lambda(3) = eigs(7);

save('lambdas', 'lambda');
% save('eigs', 'eigs');

end


