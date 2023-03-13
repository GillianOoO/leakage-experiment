function [res] = In_compute_num(num)
% if num belongs to computational subspace return 1, num in {0,...3^Nq -1}
% else return 0.

arr_3base_num = dec2base(num,3);
res = 1;
if length(arr_3base_num) == 1
    if arr_3base_num(1) == '2'
        res = 0;
        return
    else
        return
    end
end
for k = 1 : length(arr_3base_num)
    if arr_3base_num(k) == '2'
        res = 0;
        return
    end
end

end