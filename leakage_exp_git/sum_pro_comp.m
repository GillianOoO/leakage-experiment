function [pro] = sum_pro_comp(diag_pro, Nq)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
% fprintf('Nq:%d\n', Nq);
MAX = 3^Nq;
pro = 0;
for j = 0:MAX-1

    j_base3_str = dec2base(j,3);
    flag=0;
    for cur_str = 1 : length(j_base3_str)
       if j_base3_str(cur_str) == '2'
           flag = 1;
           break;
       end
    end
    if flag == 0
%         pro = pro + diag_pro(j+1,j+1);
        for k = 1 : 3^Nq
            pro = pro + diag_pro(j+1,k);
        end
%         fprintf('%d-th: pro=%f\n',j,pro);
    end
end


end

