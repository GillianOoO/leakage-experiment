function [res] = If_IthqubitL(num,I,Nq)
% if only the Ith qubit of num is in leakage subspace return 1.
% arr_full: length Nq, the i-th element in {0,1,2}, denote the i-th qubit.
% else return 0. arr_3base_num: char array.
% 

arr_3base_num = dec2base(num,3);

arr_full = zeros(1,Nq);
m = length(arr_3base_num);
res = 1;
for k = 1 : m
    arr_full(k) = arr_3base_num(m-k+1)-'0';
end

% display(arr_3base_num);
% fprintf('%d ', arr_full);

if arr_full(I) ~= 2
    res = 0;
    return;
end

for k = 1 : Nq
    if k == I
        continue;
    end
    if arr_full(k) == 2
        res = 0;
        return;
    end
end

end