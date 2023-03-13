function [f] = LRB(NumRand, NumM,pro_noise)
%return 
fprintf('rand = %d, m = %d\n', NumRand, NumM);
f = zeros(1, NumM);

p1 = pro_noise(1);
p2 = pro_noise(2);
p3 = pro_noise(3);
dep_noise = pro_noise(4);
rho_ini = zeros(3,3);
rho_ini(1,1) = 1;%|0><0|


%%pzto:ZeroToOne, potz: OneToZero, pl: leak, ps: seepage
pzto=0.05;
potz=0.05;
pl=0.0005;
ps=0.0005;
meas_noise = [1-pzto, pzto, 0;potz, 1-potz-pl, pl; 0, ps, 1-ps].';

for cur_m = 1: NumM
   f(cur_m) = 0;
   for cur_rand = 1 : NumRand
       rho = rho_ini;
       for cur_gate = 1 : cur_m
           rho = Gen_noise(rho, dep_noise, p1, p2, p3);
           gate = Gen_pauli_gate();
           rho = gate' * rho * gate; 
%            fprintf('cur_gate=%d, rho=',cur_gate);
%            display(rho);
       end
       %add measurement noise:
       diag_rho = diag(rho);
       diag_rho = meas_noise * diag_rho;
       %fidelity
       f(cur_m)= f(cur_m) + diag_rho(1)+diag_rho(2);
    %   fprintf('cur_rand=%d, f:%f\n', cur_rand, rho(1,1)+rho(2,2));
   end

   f(cur_m) = f(cur_m)/NumRand;
end






end

