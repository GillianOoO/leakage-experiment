function rho_new = Gen_noise(rho, p, p1, p2, p3)
%UNTITLED3 generate depolarize noise with pro p and leakage L
%the value for the specific elements: [1,2]:p_1,[2,1]:p_2,[2,2]:p_3
%p: depolarizing noise.
%1-qubit

leak_mat = [0 0 0;0 0 0;0 0 p3];

identity = zeros(3,3);
identity(1:2,1:2) = eye(2);
rho_new = (1-p3)*((1-p)*rho + p/2 * identity) + leak_mat;



end

