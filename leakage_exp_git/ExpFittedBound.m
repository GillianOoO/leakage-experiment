clear all;

%estimate the upper and lower bound of the fitted single exponential decay
%curve.

Y0_lower= 0.999;
Y0_upper = 0.9991;
Plateau_lower = -1.96;
Plateau_upper = 0.9762;
K_lower = 3.864e-6;
K_upper = 0.0005548;
K = 0.0002793;

Nq = 3;

file = strcat('FidBitDamp_', num2str(Nq));
file = strcat(file, '.txt');

[in_data]=load(file);

in = in_data(:,1);


%Y_lower=(Y0_lower - Plateau_upper)*exp(-K_upper*in) + Plateau_lower*ones(size(in));
%Y_upper=(Y0_upper - Plateau_lower)*exp(-K_lower*in) + Plateau_upper*ones(size(in));

lambda_lower = exp(-K_upper);

lambda_upper = exp(-K_lower);

lambda = exp(-K);


fprintf('lambda: %f, lambda_lower: %f, lambda_upper: %f\n', lambda, lambda_lower, lambda_upper);

