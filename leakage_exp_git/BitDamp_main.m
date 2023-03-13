clear all;

p_meas=[0.05,0.1];
pzto=p_meas(1);%0.05;
potz=p_meas(2);%0.1;


pl= 0.0001;%0.000008;%0.005;
ps= pl;%0.005;
pm_l0 = 0.0001;
pm_s0 = 0.0001;
pm_l1 = 0.0005;
pm_s1 = 0.0005;


meas_noise = [1-pzto - pm_l0, pzto, pm_l0;potz, 1-potz-pm_l1, pm_l1; pm_s0, pm_s1, 1-pm_s0-pm_s1].';

% step=90, NumRepeat=200, NumM=9000, step = 101, NumRepeat = 100,
% NumM=5000, eps = 0, pl=ps=0.0001
% step = 50, NumRepeat = 100, NumM=5000, eps = 0.0005, pl=ps=0
Nq = 4;
step = 50; 
NumRepeat = 200;
NumM = 2000;%20000

ErrVal = load('ErrVal.txt');%generated by 'BitDamp_gen_noise.m'
ErrVal = ErrVal(:, 1: Nq);
% EBitDamp = eye(L_mat);
%%%vec_a = [11,11,110,1100,...,110...0]
%%%vec_b = [02,20,200,2000,...,200...0]

start = tic;
[f, err_bar] = LRB_BitDamp(Nq, step, NumRepeat, NumM,pl, ps, meas_noise, ErrVal);
time = toc(start);

fprintf('Time cost: %f\n', time);
cur_f = 1;
for k = 1 : step: NumM
    fprintf('%d %f %f\n',k, f(cur_f), err_bar(cur_f));
    cur_f = cur_f + 1;
end
in = 1 :step: NumM;

out_file = strcat('FidBitDamp_', num2str(Nq));
out_file = strcat(out_file, '.txt');
fp = fopen(out_file,'w+'); 
f = f';

cur_f = 1;
for k = 1 : step: NumM
    fprintf(fp, '%d %f %f\n',k, f(cur_f), err_bar(cur_f));
    cur_f = cur_f + 1;
end

% x = in;
% g = fittype('a1+a2*exp(-a3*x)+a4*exp(-a5*x)+a6*exp(-a7x)','coefficient',{'a1','a2', 'a3', 'a4', 'a5', 'a6','a7'});
% [f0,gof,output] = fit(in.',f.',g,'StartPoint',[[ones(size(x)); exp(-x)].'\f.';1]);
% xx = linspace(1,NumM,length(x));
% 
% hold on
% p1 = plot(in, f, 'o','LineWidth',2);
% pfit=plot(xx, f0(xx),'-','LineWidth',2);
% hold off


