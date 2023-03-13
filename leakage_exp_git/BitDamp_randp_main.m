clear all;

%random generate p_ij: leak p_ij smaller, and same subspace large.

Nq = 3;
ExpN = 3^Nq;
noise_strength = [1e-3,1e-6];% [1e-4 (1 pm precision),1e-7(1 pm precision)];
precision = 5e-1;
%get pp in 3^Nq * 3^Nq, p_i in 2 * Nq, ErrorChannel(:,:,i,j) denotes Eij,
%and ErrorChannel(:,:,ExpN,ExpNq)
[pp,p, ErrorChannel] = BitDamp_randp_noiseChannelGen(Nq, ExpN,noise_strength,precision);


save('pp_randp_SPFree.txt','pp', '-ascii');
save('p_randp_SPFree.txt','p', '-ascii');

%generate circuit
p_meas=[0.05,0.1];
pzto=p_meas(1);%0.05;
potz=p_meas(2);%0.1;
spl = 0;%0.0001;
sps = 0;%0.0001;
pm_l0 = 0.0001;
pm_s0 = 0.0001;
pm_l1 = 0.0005;
pm_s1 = 0.0005;
meas_noise = [1-pzto - pm_l0, pzto, pm_l0;potz, 1-potz-pm_l1, pm_l1; pm_s0, pm_s1, 1-pm_s0-pm_s1].';

% step = 50; 
% NumRepeat = 200;
% NumM = 2000;
step = 50; 
NumRepeat = 200;
NumM = 2000;

start = tic;
[f, err_bar] = LRB_BitDamp_randp(Nq, step, NumRepeat, NumM, spl, sps, meas_noise, ErrorChannel, pp);
time = toc(start);
fprintf('Time cost: %f\n', time);

cur_f = 1;
for k = 1 : step: NumM
    fprintf('%d %f %f\n',k, f(cur_f), err_bar(cur_f));
    cur_f = cur_f + 1;
end
in = 1 :step: NumM;

out_file = strcat('FidBitDamprandp_SPFree_', num2str(Nq));
out_file = strcat(out_file, '.txt');
fp = fopen(out_file,'w+'); 
f = f';

cur_f = 1;
for k = 1 : step: NumM
    fprintf(fp, '%d %f %f\n',k, f(cur_f), err_bar(cur_f));
    cur_f = cur_f + 1;
end


