clear all;

K_pauli = 1.989e-5;
K_lower_pauli = 1.902e-5;
K_upper_pauli = 2.075e-5;

K_iswap = 0.0002185;
K_lower_iswap = 0.0002163;
K_upper_iswap = 0.0002208;


lambda_lower_pauli = exp(-K_upper_pauli);
lambda_upper_pauli = exp(-K_lower_pauli);
lambda_pauli = exp(-K_pauli);

lambda_lower_iswap = exp(-K_upper_iswap);
lambda_upper_iswap = exp(-K_lower_iswap);
lambda_iswap = exp(-K_iswap);

p_pauli = (1- lambda_pauli)/4;
p_upper_pauli = (1 - lambda_lower_pauli)/4;
p_lower_pauli = (1 - lambda_upper_pauli)/4;

% res = 1 - lambda_iswap/lambda_pauli;
% res_upper = 1 - lambda_lower_iswap/lambda_upper_pauli;

res = (lambda_pauli - lambda_iswap)/(3*lambda_pauli-2);
res_upper = ( lambda_upper_pauli- lambda_lower_iswap)/(3*lambda_upper_pauli-2);

% res = ((1 - lambda_iswap) - p_pauli)/(1 - 3 * p_pauli);
% res_upper = (1 - lambda_lower_iswap - p_lower_pauli)/(1 - 3 * p_upper_pauli);

leakage = res/2;
seepage = res*2/5;

fprintf('lambda_iswap: %fpm%f\nlambda_pauli:%fpm%f\n',lambda_iswap,...
    lambda_upper_iswap -lambda_iswap, lambda_pauli, lambda_upper_pauli-lambda_pauli);

fprintf('L=%fpm%f,S=%fpm%f\n',leakage,(res_upper - res)/2,...
    seepage,(res_upper-res)*2/5);

fprintf('Pauli:%fpm%f\n', p_pauli, p_upper_pauli -p_pauli);


