function [f, err_bar] = LRB_BitDamp_randp(Nq, step, NumRepeat, NumM, pl, ps, meas_noise, ErrorChannel, pp)
% NumRepeat: the number of repetitions
% NumM: length of the circuit
% pp:3^n by 3^n matrix pp_ij: pro from i to j.

f = zeros(1, length(1:step:NumM));

I1_single = [1,0,0;0,1,0;0,0,0];
I1 = eye(3^Nq);

for k = 1 : Nq
    I1 = I1 *kron(eye(3^(k-1)), kron(I1_single, eye(3^(Nq - k))));  
end

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


%%%% there are cur_m gate: U_1,...,U_m in total
cur_f = 1;
err_bar = zeros(size(1:step:NumM));
for cur_m = 1: step: NumM
    pro_comp = zeros(1,NumRepeat);
    for cur_rand = 1 : NumRepeat
        rho = rho_ini;
        rho = (1 - pl - ps) * rho + ps * I1/dim_I1 + pl * I2/dim_I2;%pre error
        for cur_gate = 1 : cur_m        
            % perform noise + Pauli gate 
            rho = BitDampNoise_randp(rho, ErrorChannel, pp);
            gate_multi = sparse(eye(3^Nq));
            for cur_q = 1 : Nq
                gate = Gen_pauli_gate();
                gate_multi = gate_multi * sparse(kron(eye(3^(cur_q-1)),sparse(kron(gate,eye(3^(Nq-cur_q))))));
            end
            rho = gate_multi * rho * gate_multi';
            %perform noise

%             fprintf([repmat('%f ', 1, size(rho,2)) '\n'], rho');
%             fprintf('%d, %d, %d) trace of rho:%f\n', cur_m,cur_rand, cur_gate, trace(rho));
        end
    
%         fmt=[repmat('%g+%gi\t',1,size(rho,2)-1) '%g+%gi\n'];
%         for i=1:size(rho,1)
%             fprintf(fmt,reshape(reshape([real(rho(i,:)),imag(rho(i,:))],size(rho,2),[]).',1,[]));
%         end
%         fprintf([repmat('%f ', 1, size(rho,2)) '\n'], rho');
        %add measurement noise:
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