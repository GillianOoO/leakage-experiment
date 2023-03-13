function [pauli_gate] = Gen_pauli_gate()
%UNTITLED2 generate a gate in pauli group

% sign_mat = [1,-1,1i,-1i];
% sign_ind = randi(4);


pauli(:,:,1) = [0 1;1 0];
pauli(:,:,2) = [0 -1i; 1i 0];
pauli(:,:,3) = [1 0; 0 -1];
pauli(:,:,4) = eye(2);

pauli_ind = floor(randi(4));
% pauli_gate_ini = sign_mat(sign_ind) * pauli(:,:,pauli_ind);
pauli_gate_ini = pauli(:,:,pauli_ind);

pauli_gate = eye(3);
pauli_gate(1:2,1:2)=pauli_gate_ini;

end

