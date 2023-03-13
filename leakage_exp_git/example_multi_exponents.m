clear all;
% get data
% dx = 0.02;
% x  = (dx:dx:1.5)';
% y  =  5*exp(0.5*x) + 4*exp(-3*x) + 2*exp(-2*x) - 3*exp(0.15*x);
step = 50;
eps = 0.0005;
pl = 0.0001;
ps = pl;

data_file = strcat('fidelity', num2str(step));
data_file = strcat(data_file, '_');
data_file = strcat(data_file, num2str(eps));
data_file = strcat(data_file, '_');
data_file = strcat(data_file, num2str(pl));
data_file = strcat(data_file, '.txt');

[data] = load(data_file);

x = data(:,1);
y_ori = data(:,2);

% x = 1:50:10000;
% x = x';

eps = 0.0005;
pl = 0.0001;
ps = pl;
a5 = 1-pl-ps;
a6 = a5*(1-eps/4);
a7 = a5*(1-eps/2);


a1 = 0.4740657;
a2 = 0.03631479;
a3 = 0.02256783;
a4 = 0.4667311;
y = a1+a2*a5.^x+a3*a6.^x+a4*a7.^x;
fprintf('Fit_3phase:\n%f',a1);
fprintf('+%f*%f^x', a2, a5);
fprintf('+%f*%f^x', a3, a6);
fprintf('+%f*%f^x\n', a4, a7);

y=y_ori;
Y0=0.9997;
Plateau = 0.4708;
PercentFast = 93.34;
KFast = 0.0004444;
KSlow = 0.0001491;

SpanFast=(Y0-Plateau)*PercentFast*.01;

SpanSlow=(Y0-Plateau)*(100-PercentFast)*.01;

y_twophase=Plateau + SpanFast*exp(-KFast*x) + SpanSlow*exp(-KSlow*x);
fprintf('Fit_2phase:\n%d',Plateau);
fprintf('+%d*%d^x', SpanFast, exp(-KFast));
fprintf('+%d*%d^x\n', SpanSlow, exp(-KSlow));


% calculate integrals
iy1 = cumtrapz(x, y);
iy2 = cumtrapz(x, iy1);
iy3 = cumtrapz(x, iy2);

% get exponentials lambdas
Y = [iy1, iy2, iy3, x.^3, x.^2, x, ones(size(x))];
A = pinv(Y)*y;

lambdas = eig([A(1), A(2), A(3); 1, 0, 0; 0, 1, 0]);
% lambdas


% get exponentials multipliers
X = [ones(size(x)), exp(lambdas(1)*x), exp(lambdas(2)*x), exp(lambdas(3)*x)];
P = pinv(X)*y;
yfit = P(1)+P(2)*exp(lambdas(1)*x)+P(3)*exp(lambdas(2)*x)+P(4)*exp(lambdas(3)*x);
% P
hold on
p0 = plot(x, y_ori, 'x','LineWidth',2);
p1 = plot(x, y, 'd','LineWidth',2);
p2 = plot(x, y_twophase, '--','LineWidth',2);
%pfit=plot(x, f0(x),'-','LineWidth',2);
%plinear=plot(x, yfit,'-.','LineWidth',2);
hold off
%legend('Simulated', 'Fitted', 'FittedNew');
legend('ori', 'Fitted3phase','Fitted2phase');

% exponents = zeros(3);
% fprintf('%f',P(1));
% for k = 1 : 3
%     fprintf('+%f*%f^x', P(k+1), exp(lambdas(k)));
% end