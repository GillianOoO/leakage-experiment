function [f,err_bar] = ILRB_pauli(step, NumRepeat, NumM,p_pauli, pl, ps, meas_noise)
% NumRepeat: the number of repeations
% NumM: length of the circuit

Nq = 2;
L_mat = Nq;
I1 = [1, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0;
    ];
I = eye(3^Nq);
I2 = I - I1;

dim_I1 = 2^Nq;
dim_I2 = 3^Nq - 2^Nq;

phi = zeros(3^Nq,1);
phi(1) = 1;
rho_ini = sparse(phi * phi');

meas_noise_multi = eye(3^Nq);
for j = 1 : Nq
    meas_noise_multi = sparse(meas_noise_multi) * sparse(kron(eye(3^(j-1)),sparse(kron(meas_noise, eye(3^(Nq-j))))));
end

%%%vec_a = [11,11]
%%%vec_b = [02,20]
vec_a = zeros(L_mat,1);%column vector, represent the flip bit.
vec_b = zeros(L_mat,1);
Err_pauli = p_pauli * ones(2,Nq);
%%%% there are cur_m gate: T,U1,T,...,U_m in total
for k = 1 : Nq
    if k == 1
        ind_a = 3^k + 1 + 1;
    else
        ind_a = 3^(k-2) + 3^(k-1)+1;
    end
    ind_b = 2*3^(k-1) + 1;
    vec_a(k) = ind_a;
    vec_b(k) = ind_b;
end

f = zeros(1, length(1:step:NumM));
err_bar = zeros(size(1:step:NumM));


cur_f = 1;
for cur_m = 1: step: NumM
    pro_comp = zeros(1,NumRepeat);
    for cur_rand = 1 : NumRepeat
        rho = rho_ini;
        rho = (1 - pl - ps) * rho + ps * I1/dim_I1 + pl * I2/dim_I2;%pre error
        for cur_gate = 1 : cur_m  
            %perform noise for pauli group
            rho = BitDampNoise(rho, Err_pauli, vec_a, vec_b);
%             for cur_q = 1 : Nq
%                 rho = rho -  p_pauli * rho(vec_a(cur_q), vec_a(cur_q)) * vec_rep(vec_a(cur_q),L_mat) * vec_rep(vec_a(cur_q),L_mat)'...
%                     - p_pauli * rho(vec_b(cur_q), vec_b(cur_q)) * vec_rep(vec_b(cur_q),L_mat) * vec_rep(vec_b(cur_q),L_mat)'...
%                    + p_pauli * rho(vec_a(cur_q), vec_a(cur_q)) * vec_rep(vec_b(cur_q),L_mat) * vec_rep(vec_b(cur_q),L_mat)'...
%                      + p_pauli * rho(vec_b(cur_q), vec_b(cur_q)) * vec_rep(vec_a(cur_q),L_mat) * vec_rep(vec_a(cur_q),L_mat)';
%             end    
            gate_multi = eye(3^Nq);
            for cur_q = 1 : Nq
                gate = Gen_pauli_gate();
                gate_multi = gate_multi * sparse(kron(eye(3^(cur_q-1)),sparse(kron(gate,eye(3^(Nq-cur_q))))));
            end
            rho = gate_multi * rho * gate_multi';
        end
        diag_rho = diag(rho);%a 2^Nq length vector
        diag_rho = diag(diag_rho);% a 2^Nq by 2^Nq matrix
        diag_rho = meas_noise_multi * diag_rho;
        %fidelity
        pro_comp(cur_rand) = sum_pro_comp(diag_rho, Nq);
        f(cur_f)= f(cur_f) + pro_comp(cur_rand);
        %   fprintf('cur_rand=%d, f:%f\n', cur_rand, rho(1,1)+rho(2,2));
    end
    f(cur_f) = f(cur_f)/NumRepeat;
    err_bar(cur_f) = BitDamp_errbar_cal(f(cur_f), pro_comp);
    cur_f = cur_f + 1;
end


end