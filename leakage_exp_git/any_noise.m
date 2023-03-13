clear all;

p_meas=[0.05,0.1];
pzto=p_meas(1);%0.05;
potz=p_meas(2);%0.1;

eps = 0.00028;%0.000082;%0.00028;
pl= 0;%0.000008;%0.005;
ps= pl;%0.005;
pm_l = 0.0028;
pm_s = 0.0028;

meas_noise = [1-pzto, pzto, 0;potz, 1-potz-pm_l, pm_l; 0, pm_s, 1-pm_s].';


step = 50; 
% step=90, NumRepeat=200, NumM=9000, step = 101, NumRepeat = 100,
% NumM=5000, eps = 0, pl=ps=0.0001
% step = 50, NumRepeat = 100, NumM=5000, eps = 0.0005, pl=ps=0

NumRepeat = 40;
NumM = 1000;%20000
start = tic;
[f] = LRB_any_noise(step, NumRepeat, NumM, pl, ps, meas_noise);
time = toc(start);

fprintf('Time cost: %f\n', time);

in = 1 :step: NumM;
f = abs(f);
plot(in,f,'-o','LineWidth',2);
xlabel('Depth');
ylabel('Pr');
set(gca,'Fontsize',18);
out_file = strcat('FidelityGenNoise', num2str(step));
out_file = strcat(out_file, '.txt');
fp = fopen(out_file,'w+'); 
f = f';

cur_f = 1;
for k = 1 : step: NumM
    fprintf(fp, '%d %f\n',k, f(cur_f));
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


