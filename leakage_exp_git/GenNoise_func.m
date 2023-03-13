function [rho_new] = GenNoise_func(rho)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

G(:,:,1) = eye(3);
G(:,:,2) = [1 0 0;
            0 0 1;
            0 1 0];
G(:,:,3) = [0 0 1;
            0 1 0;
            1 0 0];
G(:,:,4) = [1 0 0;
            0 exp(2*pi*1i/3) 0;
            0 0 exp(4*pi*1i/3)];
G(:,:,5) = [1 0 0;
            0 exp(4*pi*1i/3) 0;
            0 0 exp(2*pi*1i/3)];

M = 5;

rho_new = zeros(size(rho));
p_noise = rand(1,M)/100;
p_noise(1) = 1 - sum(p_noise(2:M));
for k = 1 : M
    rho_new = rho_new + p_noise(k) * kron(G(:,:,k),G(:,:,1)) * rho * kron(G(:,:,k),G(:,:,1))';
end

p_noise = rand(1,M)/100;
p_noise(1) = 1 - sum(p_noise(2:M));
rho_new = zeros(size(rho));
for k = 1 : M
    rho_new =  rho_new + p_noise(k) * kron(G(:,:,1),G(:,:,k)) * rho * kron(G(:,:,1),G(:,:,k))';
end

        
end

