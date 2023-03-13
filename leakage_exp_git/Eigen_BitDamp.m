clear all;

ErrVal = load('ErrVal.txt');

Nq = size(ErrVal,2);
mat = zeros(Nq+1, Nq+1);

sum_Halfp = 0;
for k = 1 : Nq
    mat(k+1,1) = ErrVal(1,k)/2;
    sum_Halfp = sum_Halfp + ErrVal(1,k)/2;
    mat(1,k+1) = ErrVal(2,k);
    mat(k+1,k+1)= 1- ErrVal(2,k);
end

mat(1,1) = 1 - sum_Halfp;

% for k = 1 : Nq+1
%     fprintf('%e ', mat(k,:));
%     fprintf('\n');
% end

eigs = eig(mat);

fprintf('%e ', eigs);
