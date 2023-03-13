function [f] = LRB_any_noise(step, NumRepeat, NumM, pl, ps, meas_noise)
% NumRepeat: the number of repeations
% NumM: length of the circuit

Nq = 2;
f = zeros(1, length(1:step:NumM));

% Enoise = [
%     1, 0, 0,        0, 0,       0, 0,       0, 0;
%     0, 0, 0,        1, 0,       0, 0,       0, 0;
%     0, 0, 1-eps/2,  0, eps/2,   0, 0,       0, 0;
%     0, 1, 0,        0, 0,       0, 0,       0, 0;
%     0, 0, eps/2,    0, 1-eps,   0, eps/2,   0, 0;
%     0, 0, 0,        0, 0,       1, 0,       0, 0;
%     0, 0, 0,        0, eps/2,   0, 1-eps/2, 0, 0;
%     0, 0, 0,        0, 0,       0, 0,       1, 0;
%     0, 0, 0,        0, 0,       0, 0,       0, 1
%     ];


Dim = 9;
Enoise = zeros(Dim, Dim);
for k = 1 : Dim
    Enoise_cur = rand(1,Dim)/10000;
    Enoise(k,:) = Enoise_cur;
    sum_row = sum(Enoise_cur)-Enoise_cur(1,k);
    Enoise(k,k) = 1 - sum_row;
end


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

phi = zeros(3^Nq,1);
phi(1) = 1;
rho_ini = (phi * phi');

meas_noise_multi = eye(3^Nq);
for j = 1 : Nq
    meas_noise_multi = (meas_noise_multi) * (kron(eye(3^(j-1)),(kron(meas_noise, eye(3^(Nq-j))))));
end


%%%% there are cur_m gate: (U_0U_1),...,U_m in total, and U_k = P * Ed *
%%%% Enoise

cur_f = 1;
for cur_m = 1: step: NumM
%     f(cur_f) = 0;
    for cur_rand = 1 : NumRepeat
        rho = rho_ini;
 %       gate_last_multi = (eye(3^Nq));
        
        for cur_gate = 1 : cur_m
            %            if cur_gate > 1 %cur_gate = 1: combine U_0 and U_1.
            %                 rho = (Gen_noise_x(rho, p_x, Nq));
            %            end
%             rho = diag(rho);
%             rho = (Enoise * rho);
%             rho = diag(rho);
            rho = GenNoise_func(rho);
            rho = (1 - pl - ps) * rho + ps * I1/4 + pl * I2/5;
            gate_multi = (eye(3^Nq));
%             fmt=[repmat('%g+%gi\t',1,size(rho,2)-1) '%g+%gi\n'];
%             for i=1:size(rho,1)
%                 fprintf(fmt,reshape(reshape([real(rho(i,:)),imag(rho(i,:))],size(rho,2),[]).',1,[]));
%             end
            for cur_q = 1 : Nq
                gate = Gen_pauli_gate();
                gate_multi = gate_multi * (kron(eye(3^(cur_q-1)),(kron(gate,eye(3^(Nq-cur_q))))));
            end
     %       gate_last_multi =  gate_multi * gate_last_multi;
            rho = gate_multi * rho * gate_multi';
            rho = GenNoise_func(rho);
        end
        
 %       rho = (1 - pl - ps) * rho + ps * I1/4 + pl * I2/5;
 %       rho = gate_last_multi' * rho * gate_last_multi;
%         fmt=[repmat('%g+%gi\t',1,size(gate_last_multi,2)-1) '%g+%gi\n'];
%         for i=1:size(gate_last_multi,1)
%             fprintf(fmt,reshape(reshape([real(gate_last_multi(i,:)),imag(gate_last_multi(i,:))],size(gate_last_multi,2),[]).',1,[]));
%         end
%         fmt=[repmat('%g+%gi\t',1,size(rho,2)-1) '%g+%gi\n'];
%         for i=1:size(rho,1)
%             fprintf(fmt,reshape(reshape([real(rho(i,:)),imag(rho(i,:))],size(rho,2),[]).',1,[]));
%         end
%         fprintf([repmat('%f ', 1, size(rho,2)) '\n'], rho');
        %add measurement noise:
        % fprintf('(m=%d, rand=%d)---\n',cur_m, cur_rand);
        diag_rho = diag(rho);%a 2^Nq length vector
        diag_rho = diag(diag_rho);% a 2^Nq by 2^Nq matrix
        %        display(rho);
        %        display(diag_rho);
        %        display(meas_noise_multi);
        diag_rho = meas_noise_multi * diag_rho;
        %fidelity
%         display(diag_rho);
        pro_comp = sum_pro_comp(diag_rho, Nq);
        f(cur_f)= f(cur_f) + pro_comp;
        %   fprintf('cur_rand=%d, f:%f\n', cur_rand, rho(1,1)+rho(2,2));
    end
    
    f(cur_f) = f(cur_f)/NumRepeat;
    cur_f = cur_f + 1;
end


end