% generate pauli group: Gen_pauli_gate()
% generate error model: rho_new = Gen_noise(rho, p, L)

clear all;

% rho = eye(3)/2;
% rho(3,3) = 0;
% 
% pauli = Gen_pauli_gate();
% display(pauli);
% 

% 
% rho_new = Gen_noise(rho, p, p1, p2, p3);
% 
% display(rho_new);

% dep_noise = 0.02;
% p1=0.001;
% p2=0.001;
% p3=0.001;
start = tic;

NumRand = 30;
NumM=40;


pro_noise = [0.005,0.005]; % [leak,leak,leak,dep_noise]

l1 = pro_noise(1);
%l2 = 1-(1-pro_noise(3))*(1-pro_noise(4))-pro_noise(3);
l2=pro_noise(2);

%%pzto:ZeroToOne, potz: OneToZero, pl: leak, ps: seepage
p_meas=[0.05,0.1,0.05,0.05];

pzto=p_meas(1);%0.05;
potz=p_meas(2);%0.1;
pl=p_meas(3);%0.005;
ps=p_meas(4);%0.005;
% meas_noise = [1-pzto, pzto, 0;potz, 1-potz-pl, pl; 0, ps, 1-ps].';
meas_noise = [1-pzto, pzto, 0;potz, 1-potz-pl, pl; 0, ps, 1-ps].';

in = 1:NumM;
Nq =2;

lm = (2^Nq - (2 - pl)^Nq)/2^Nq;
sm = ((2-pl+ps)^Nq - (2 - pl)^Nq)/(3^Nq - 2^Nq);

a = (l2-l2*lm+l1*sm)/(l1+l2);
% a = l2/(l1+l2);
b=l1*(1-lm-sm)/(l1+l2);

lambda=1-l1-l2;

f = LRB_multi(NumRand, NumM, pro_noise, Nq, meas_noise);
    
% display(f);

cal_cur = a * ones(1,NumM) + b *(lambda).^(in+ones(1,length(in)));
% cal_cur = a * ones(1,NumM) + b * lambda.^in;
x = in;
g = fittype('a+b*exp(-c*x)');
f0 = fit(in.',f.',g,'StartPoint',[[ones(size(x)); exp(-x)].'\f.';1]);
xx = linspace(1,NumM,length(x));

hold on
p1 = plot(in, f, 'o','LineWidth',2);
pfit=plot(xx, f0(xx),'-','LineWidth',2);
p2=plot(in, cal_cur,'-.*', 'LineWidth',1);
hold off


% p1.marker('sq');
% p2.Marker='sq';
xlabel('The number of gates','fontsize',18);
ylabel('Fidelity','fontsize',18);
legend('Simulated','fitted curve', 'Theoretical','fontsize',18);
title('')




a_fit = f0.a;
lam_fit = exp(-1 * f0.c);
b_fit = f0.b/lam_fit;

leak_fit_a = b_fit * (1 - lam_fit)/(1 - lm - sm);
seep_fit_a = 1 - lam_fit - leak_fit_a;

leak_fit_b = (1-lam_fit)*(1- lm - a)/(1-lm-sm);
seep_fit_b = 1 - lam_fit - leak_fit_b;

% leak_fit_m = b_fit * (1-lm) * (1-lam_fit)/(b_fit + a_fit)/(1-lm-sm);
% seep_fit_m = 1 - lam_fit - leak_fit_m;

time = toc(start);
fprintf('Qubit=%d, time=%.2f\n', Nq, time);
fprintf('theoretical: a = %.4f, b = %.4f, lam=%.4f\n', a, b, lambda);
fprintf('Fitted: a = %.4f, b = %.4f, lam=%.4f\n', a_fit, b_fit, lam_fit);
fprintf('L1=%.4f, L2=%.4f.\nFitted: leak_a=%.4f, seep_a=%.4f, leak_b=%.4f, seep_b=%.4f\n',...
    l1,l2, leak_fit_a, seep_fit_a, leak_fit_b, seep_fit_b);
% display(f0);
% fprintf('leak_c = %.4f, seep_c = %.4f\n', leak_fit_m, seep_fit_m);




