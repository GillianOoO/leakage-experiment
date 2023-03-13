clear all;

Nq = 3;

file = strcat('FidBitDamprandp_', num2str(Nq));
file = strcat(file, '.txt');

[in_data]=load(file);

in = in_data(:,1);
f = in_data(:,2);
err = in_data(:,3);

%fit single exponential decay
%fit with prism 9.

%%%for four qubits.
% Y0 = 0.9987;
% Plateau = 0.2164;
% K = 4.286e-5;

%%%for three qubits randp exp.
Y0 = 0.999;
Plateau = 0.4190;
K = 2.597e-5;

Y=(Y0 - Plateau)*exp(-K*in) + Plateau*ones(size(in));
figure
box on
hold on
fill([in;flipud(in)],[f-err;flipud(f+err)],[0.70 0.80 0.940],'linestyle','none');
s1=line(in,f,'LineStyle',':','Color',[0 0.4470 0.7410],'marker','d','linewidth',2);

s2=plot(in,Y,'-','Color',[0.4660 0.6740 0.1880],'linewidth',2);
legend('','Experiment','Fitted','fontsize',18);

xlabel('Size','fontsize',18);
ylabel('Pr','fontsize',18);
set(gca,'fontsize',12);