function [rho_new] = BitDampNoise(rho, ErrVal, vec_a, vec_b)
%performing the leakage bit damp noise on rho, and return updated rho
% rho_new = E_i rho E_i^dagger
%ErrVal(1,i): p_i
%ErrVal(2,i): q_i
%i: 11->02, 11->20, 110->200, ...

%E_{2M + 1} = ...
%p(j) = p(vec_a(j));
%E_{vec_a(i),vec_b(i)} = sqrt(ErrVal_..)ket(vec_b(i)) bra(vec_a(i))

M = length(ErrVal); %2 * M: the number of damps. M=Nq

ExpNq = length(rho);

Error = zeros(ExpNq, ExpNq, 2 * M + 1);

Error(:,:,2*M + 1) = eye(ExpNq);

%a_k = {11,11,110,1100,11000,...}
%b_k = {02,20,220,2200,22000,...}

for k = 1 : M
   Error(vec_a(k),vec_a(k), 2 * M + 1) = sqrt(1 - ErrVal(1,k));  
   Error(vec_b(k),vec_b(k), 2 * M + 1) = sqrt(1 - ErrVal(2,k));
   Error(vec_b(k), vec_a(k), 2 * k - 1) = sqrt(ErrVal(1,k));
   Error(vec_a(k), vec_b(k), 2 * k) = sqrt(ErrVal(2,k));
end

Error(vec_a(1),vec_a(1), 2 * M + 1) = sqrt(1 - ErrVal(1,1) - ErrVal(1,2));
%note that a(1) = a(2)

rho_new = zeros(size(rho));

for k = 1 : 2 * M + 1
   rho_new = rho_new + Error(:,:,k) * rho * Error(:,:,k)';
end

end

