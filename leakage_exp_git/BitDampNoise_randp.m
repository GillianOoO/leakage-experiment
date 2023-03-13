function [rho_new] = BitDampNoise_randp(rho, ErrorChannel, pp)
%performing the leakage bit damp noise on rho, and return updated rho
% rho_new = E_ij rho E_ij^dagger

ExpN = length(rho);
rho_new = zeros(size(rho));


for k = 1 : ExpN
    for j = 1 : ExpN
        if pp(k,j) ~= 0
            rho_new = rho_new + ErrorChannel(:,:,k,j) * rho * ErrorChannel(:,:,k,j)';
        end
    end
end

rho_new = rho_new + ErrorChannel(:,:,ExpN,ExpN) * rho * ErrorChannel(:,:,ExpN,ExpN)';

end

