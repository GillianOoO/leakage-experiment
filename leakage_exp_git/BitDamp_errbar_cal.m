function [err_bar] = BitDamp_errbar_cal(mean, pro_arr)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

len = length(pro_arr);

err_bar = 0;
for k = 1 : len
    err_bar = err_bar + (pro_arr(k) - mean)^2;
end
err_bar = sqrt(err_bar);

err_bar = 3 * err_bar/len;


end

