clear all;

Nq = 3;
K = 2.597e-5;
K_lower = 2.552e-5;
K_upper = 3.643e-5;



% Nq = 4;
% K = 4.286e-5;
% K_lower = 3.087e-5;
% K_upper = 5.484e-5;


% Nq = 5;
% K = 2.103e-5;
% K_lower = 1.312e-5;
% K_upper = 2.895e-5;



lambda = exp(-K);
lambda_lower = exp(-K_upper);
lambda_upper = exp(-K_lower);

err_1 = (lambda_upper - lambda);
err_2 = (lambda - lambda_lower);
err_lambda = max(err_1,err_2);


ave_p = (1-lambda)/(Nq+2);

%ave_p = (1-lambda)/2; %coef=0.

ave_p_upper = (1-lambda_lower)/(Nq+2);
ave_p_lower = (1-lambda_upper)/(Nq+2);
err_p = max(ave_p_upper - ave_p, ave_p - ave_p_lower);

L_ave = Nq * ave_p;
L_ave_upper = Nq * ave_p_upper;
L_ave_lower = Nq * ave_p_lower;
err_L = max(L_ave_upper - L_ave, L_ave - L_ave_lower);

S_ave = Nq * 2^Nq * ave_p/(3^Nq - 2^Nq);
S_ave_upper = Nq * 2^Nq * ave_p_upper/(3^Nq - 2^(Nq));
S_ave_lower = Nq * 2^Nq * ave_p_lower/(3^Nq - 2^(Nq));
err_S = max(S_ave_upper - S_ave, S_ave - S_ave_lower);

fprintf('lambda: %epm %e\n ave_p:%epm%e\n l_ave: %epm%e\n s_ave:%epm%e\n',...
    lambda, err_lambda, ave_p, err_p, L_ave, err_L, S_ave, err_S);

