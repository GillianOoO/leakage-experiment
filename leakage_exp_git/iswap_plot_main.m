clear all;

step = 100;
% eps = 0.0005;
% pl = 0.0001;
% ps = pl;
eps = 0.00028;
pl= 0;%0.000008;%0.005;
ps= pl;%0.005;

% step = 2;
% eps = 0.1;
% pl = 0.01;
% ps = pl;

data_file = strcat('fidelity', num2str(step));
data_file = strcat(data_file, '_');
data_file = strcat(data_file, num2str(eps));
data_file = strcat(data_file, '_');
data_file = strcat(data_file, num2str(pl));
data_file = strcat(data_file, '.txt');
if ~isfile(data_file)
    fprintf('%s', data_file);
end
[data] = load(data_file);

x = data(:,1);
f = data(:,2);
x = x';
f = f';

lambda1 = 1-pl-ps;
 lambda2 = lambda1*(1-eps/4);
 lambda3 = lambda1*(1-eps/2);

fprintf('lambda1 = %f, lambda3 = %f\n', lambda1, lambda3);
g = fittype('a1+a2*exp(-a3*x)','coefficients',{'a1', 'a2','a3'});
startpoint = [[ones(size(x)); exp(-x)].'\f.';1];
% startpoint = [ones(size(x'));ones(size(x'))].'\f;
% 
% startpoint = zeros(2,1);
% startpoint(3) = 0;
% startpoint(1:2) = inipoint(1:2);
% % startpoint(4) = inipoint(3);

x = x';
f = f';
[f0,gof,output] = fit(x,f,g,'StartPoint',startpoint);


% g_linear = fittype('a1+a2*0.9875.^x+a3*0.9750.^x','coefficients',{'a1', 'a2','a3'});
% sp_linear = [ones(size(x')); 0.9875.^x'; 0.9750.^x'].'\f;
% [f_linear,gof_linear,output_linear] = fit(x,f,g_linear,'StartPoint',sp_linear);

% g_linear = fittype('a1+a2*0.999944.^x+a4*0.999944.^x','coefficients',{'a1', 'a2','a4'});
% % g_linear = fittype('a1+a2*0.98.^x+a3*0.9555.^x+a4*0.9310.^x','coefficients',{'a1', 'a2','a3','a4'});
% % 
% sp_linear = [ones(size(x')); lambda1.^x';lambda3.^x'].'\f;
% [f_linear,gof_linear,output_linear] = fit(x,f,g_linear,'StartPoint',sp_linear);


hold on
% p1 = plot(x, f, 'o','LineWidth',2);
pfit=plot(x, f0(x),'x','LineWidth',2);
%plinear=plot(x, f_linear(x),'d','LineWidth',2);
hold off

%h = legend('Simulated', 'Fitted', 'FittedNew');
h = legend('Simulated data', 'Fitted curve');
ylabel('Pr');
xlabel('Depth');
%set(h,'FontSize',12);

fprintf('Fitted:\n%d',f0.a1);

lam = exp(-f0.a3);
fprintf('+%d*%d.^x', f0.a2, lam);
% fprintf('+%d*%d.^x', f0.a3, exp(f0.a6));
% fprintf('+%d*%d.^x\n', f0.a4, exp(f0.a7));


% fprintf('IPartial fit:\n%d',f_linear.a1);
% 
% fprintf('+%d*%d.^x', f_linear.a2, lambda1);
% % fprintf('+%d*%d.^x', f_linear.a3, lambda2);
% fprintf('+%d*%d.^x\n', f_linear.a4, lambda3);
