% generate pauli group: Gen_pauli_gate()
% generate error model: rho_new = Gen_noise(rho, p, L)

clear all;

start = tic;

NumRand = 200;
NumM = 60;
Nq = 2;

p_x = zeros(1,Nq);
for j = 1 : Nq
    p_x(j) = 0.005; 
end

%%pzto:ZeroToOne, potz: OneToZero, pl: leak, ps: seepage
p_meas=[0.05,0.1,0.05,0.05];
pzto=p_meas(1);%0.05;
potz=p_meas(2);%0.1;
pl=p_meas(3);%0.005;
ps=p_meas(4);%0.005;
meas_noise = [1-pzto, pzto, 0;potz, 1-potz-pl, pl; 0, ps, 1-ps].';




% lm = (2^Nq - (2 - pl)^Nq)/2^Nq;
% sm = ((2-pl+ps)^Nq - (2 - pl)^Nq)/(3^Nq - 2^Nq);
% a = (l2-l2*lm+l1*sm)/(l1+l2);
% % a = l2/(l1+l2);
% b=l1*(1-lm-sm)/(l1+l2);
%lambda=1-l1-l2;

fit_lam = zeros(2,Nq);
fit_px = zeros(2, Nq);


box on
for j = 1 : Nq
    f = LRB_multi_xnoise(NumRand, NumM, p_x(j), 1, meas_noise); %independent x-noise
    in = 1:NumM;
%    cal_cur = a * ones(1,NumM) + b *(lambda).^(in+ones(1,length(in)));
    x = in;
    g = fittype('a+b*exp(-c*x)');
    [f0,gof,output] = fit(in.',f.',g,'StartPoint',[[ones(size(x)); exp(-x)].'\f.';1]);
    xx = linspace(1,NumM,length(x));
    
    hold on
    p1 = plot(in, f, 'o','LineWidth',2);
    pfit=plot(xx, f0(xx),'-x','LineWidth',2);
    xlabel('The number of gates','fontsize',18);
    ylabel('Fidelity','fontsize',18);
    legend('Experimental results','fitted curve', 'fontsize',18);
    
    coef_data = confint(f0,0.95);
    
    fit_lam(1,j) = exp(-1 * coef_data(1, 3));
    fit_lam(2,j) = exp(-1 * coef_data(2, 3));
    
    fit_px(1,j) = 2 * (1 - fit_lam(1,j))/3;
    fit_px(2,j) = 2 * (1 - fit_lam(2,j))/3; 
end


% L = 1/2^Nq * ones(1,2);
% 
% S = 1/(3^Nq - 2^Nq) * ones(1,2);

tempFirst = ones(1,2);
idfirst = 1;
tempSecond = ones(1,2);
idsecond = 1;
for j = 1 : Nq
    for cur = 1 : 2
       tempFirst(cur) = tempFirst(cur) *  (2 + 2 * fit_px(cur,j));
       tempSecond(cur) = tempSecond(cur) * (2 + fit_px(cur,j));
    end
    idfirst = idfirst * (2 + 2 * p_x(j));
    idsecond = idsecond * (2 + p_x(j));
end

L = (tempFirst - tempSecond)/2^Nq;

S = (tempFirst - tempSecond)/(3^Nq - 2^Nq) ;

l1 = (idfirst - idsecond)/2^Nq;
l2 = (idfirst - idsecond)/(3^Nq - 2^Nq);

time = toc(start);
fprintf('Qubit=%d, time=%.2f\n', Nq, time);
% fprintf('theoretical: a = %.4f, b = %.4f, lam=%.4f\n', a, b, lambda);
% fprintf('Fitted: a = %.4f, b = %.4f, lam=%.4f\n', a_fit, b_fit, lam_fit);
fprintf('L=%.4f, S=%.4f.\nFitted: leak_fitted=%.4f (%.4f,%.4f), seep_fitted =%.4f (%.4f, %.4f)\n',...
    l1,l2, (L(1)+L(2))/2, L(1), L(2), (S(1)+S(2))/2, S(1), S(2));


