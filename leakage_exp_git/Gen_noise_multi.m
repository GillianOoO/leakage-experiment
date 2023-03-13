function rho_new = Gen_noise_multi(rho, dep_p, pl,Nq)
%UNTITLED3 generate depolarize noise with pro p and leakage L
%the value for the specific elements: [1,2]:p_1,[2,1]:p_2,[2,2]:p_3
%p: depolarizing noise.
%1-qubit

% leak_mat_single = [0 0 0;0 0 0;0 0 1];
% leak_mat = eye(3^Nq);
% for j = 1 : Nq
%     leak_mat = leak_mat * kron(eye(3^(j-1)), kron( leak_mat_single, eye(3^(Nq-j)) ));
% end

identity = eye(3^Nq);
id_com = [1 0 0; 0 1 0;0 0 0];
for j = 1 : Nq
    identity = identity * kron(eye(3^(j-1)), kron(id_com , eye(3^(Nq-j)) ));
end

leak_mat = eye(3^Nq) - identity;

%rho_new = (1-pl)*((1-dep_p)*rho + dep_p/(2^Nq) * identity) + pl/Nq*leak_mat;

% rho_new = (1-dep_p-pl)*rho + dep_p/(2^Nq) * identity + pl/Nq*leak_mat;
rho_new = (1-dep_p-pl)*rho + dep_p/(2^Nq) * identity + pl/(3^Nq - 2^Nq)*leak_mat;

end

