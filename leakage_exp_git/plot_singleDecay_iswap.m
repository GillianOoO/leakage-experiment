clear all;
%generate choice_startpoint.pdf

% x = 10000:100:50000;
% eps = 0.000281317;
% 
% f = (0.5+1/(2-eps)).^-x;
% 
% plot(x,f,'-','LineWidth',2);
% legend('((1-\epsilon/2)/(1-\epsilon/4))^x');
% title('Choice of the start point');

% step = 100;
% eps = 0.0005;
% pl = 0.0001;
% ps = pl;
% 
% data_file = strcat('fidelity', num2str(step));
% data_file = strcat(data_file, '_');
% data_file = strcat(data_file, num2str(eps));
% data_file = strcat(data_file, '_');
% data_file = strcat(data_file, num2str(pl));
% data_file = strcat(data_file, '.txt');
% 
% [data] = load(data_file);
% 
% x = data(:,1);
% f = data(:,2);
% 
% plot(x,f,'-o','LineWidth',2);
% ylabel('Pr');
% xlabel('Depth');
% title('Simulation curve');


step = 50;
% eps = 0.0005;
% pl = 0.0001;
% ps = pl;
eps = 0.00028;%8.2e-05;%8.2e-05;%0.00028;
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


% Y0 = 0.9967;%0.9964;%0.9960-0.9967
% Plateau = 0.4897;%0.4893;%0.489-0.4897
% K = 0.0001518;%0.0001515;%0.0001511-0.0001518

Y0 = 0.9971;%0.9973;%0.997; %"0.9966 to 0.9974"%0.9991;%%0.9989 ~ 0.9992
Plateau = 0.4977;%0.5008;%0.4998;%"0.4994 to 0.5002";%0.4859;%%0.4850 ~ 0.4869
K = 1.391e-4;%0.0001398; %"0.0001395 to 0.0001402";%5.537e-005;%5.520e-005 ~ 5.554e-005


b = Y0 - Plateau;
a = Plateau;
lam = exp(-K);
str = sprintf('%.6f + %.6f * %.6f^x', a, b, lam);

Y=(Y0 - Plateau)*exp(-K*x) + Plateau;

b1 = 1/71*(2 - pl + 2*ps)*(2 - pl) - 1/9*ps^2;
b3 = 1/8*(2 - pl - 2*ps)*(2 - pl);

lam1 = 1-pl-ps;
lam3 = (1-pl-ps)*(1-eps/2);
f_cut = f - b1 * lam1.^x - b3 * lam3.^x;

a1 = [ones(size(x'))].'\f_cut;

ideal_f = a1 + b1 * lam1.^x + b3 * lam3.^x;


ideal_str = sprintf('%.6f + %.6f * %.6f^x', a1+b1, b3, lam3);

%ideal_str = sprintf('%f + %f * %f^x + %f * %f^x', a1, b1, lam1, b3, lam3);
figure
box on
hold on
s(1) = plot(x,f,'x--','LineWidth',2);
s(2) = plot(x,Y,'-','LineWidth',2);
% s(3) = plot(x, ideal_f, ':','LineWidth',2);
xlabel('Depth');
ylabel('Pr');
set(gca,'Fontsize',18);
legend(s, {'Experimental result',str},'Fontsize',18,'Box','off');
%legend(s, {'Experimental result',str, ideal_str},'Fontsize',12, 'Location','SouthOutside','Box','off');
%legend(s(2), str,'Fontsize',12, 'Location','SouthOutside','Box','off');
%legend(s(3), ideal_str, 'Fontsize',12, 'Location','SouthOutside','Box','off');
% ylim([0.5+0.1 ,1]);
s(1).MarkerSize = 4;

fprintf('%s\n%s\n',str, ideal_str);

hold off


