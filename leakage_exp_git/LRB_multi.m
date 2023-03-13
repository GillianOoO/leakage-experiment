function [f] = LRB_multi(NumRand, NumM,pro_noise,Nq, meas_noise)
%n-qubits lrb
% fprintf('rand = %d, m = %d\n', NumRand, NumM);
f = zeros(1, NumM);


p3 = pro_noise(1);
dep_noise = pro_noise(2);
rho_ini = zeros(3^Nq,3^Nq);
rho_ini(1,1) = 1;%|0><0|



meas_noise_multi = eye(3^Nq);
for j = 1 : Nq
   meas_noise_multi = meas_noise_multi * kron(eye(3^(j-1)),kron(meas_noise, eye(3^(Nq-j))));
end

step = 1;

%%%% there are cur_m + 1 gate: (U_0U_1),...,U_m, (U_m^dagU_m-1^dag...U_1^dag) in total.
for cur_m = 1: step: NumM
   f(cur_m) = 0;
   for cur_rand = 1 : NumRand
       rho = rho_ini;
       gate_last_multi = eye(3^Nq);
       %Lambda * U_0
       rho = Gen_noise_multi(rho, dep_noise, p3, Nq);
       gate_multi = eye(3^Nq);
       for cur_q = 1 : Nq
                gate = Gen_pauli_gate();
                gate_multi = gate_multi * kron(eye(3^(cur_q-1)),kron(gate,eye(3^(Nq-cur_q))));
       end
       rho = gate_multi * rho * gate_multi';
       for cur_gate = 1 : cur_m
           if cur_gate > 1 %cur_gate = 1: combine U_0 and U_1.
                rho = Gen_noise_multi(rho, dep_noise, p3, Nq);
           end
           gate_multi = eye(3^Nq);
           for cur_q = 1 : Nq
                gate = Gen_pauli_gate();
                gate_multi = gate_multi * kron(eye(3^(cur_q-1)),kron(gate,eye(3^(Nq-cur_q))));
           end      
           gate_last_multi = gate_last_multi * gate_multi; % to obtain U_m+1
           rho = gate_multi * rho * gate_multi';
%             fprintf('cur_gate=%d, rho=',cur_gate);
%             display(rho);
       end
       %the laste gate:
       rho = Gen_noise_multi(rho, dep_noise, p3, Nq);
       gate_last_multi = inv(gate_last_multi);
       rho = gate_last_multi * rho * gate_last_multi';
       %add measurement noise:
%        fprintf('%d-th---',cur_rand);
       diag_rho = diag(rho);%a 2^Nq length vector
       diag_rho = diag(diag_rho);% a 2^Nq by 2^Nq matrix
%        display(rho);
%        display(diag_rho);
%        display(meas_noise_multi);
       diag_rho = meas_noise_multi * diag_rho;
       %fidelity
%        display(diag_rho);
       pro_comp = sum_pro_comp(diag_rho, Nq);
       f(cur_m)= f(cur_m) + pro_comp;
    %   fprintf('cur_rand=%d, f:%f\n', cur_rand, rho(1,1)+rho(2,2));
   end

   f(cur_m) = f(cur_m)/NumRand;
end






end

