clear all;


file_iswap = 'Fidelity_ILRB_iswap.txt';
file_pauli = 'Fidelity_ILRB_pauli.txt';


iswap_data=load(file_iswap);
pauli_data=load(file_pauli);

iswap_in = iswap_data(:,1);
iswap_f = iswap_data(:,2);
iswap_err = iswap_data(:,3);

pauli_in = pauli_data(:,1);
pauli_f = pauli_data(:,2);
pauli_err = pauli_data(:,3);

%fit single exponential decay
%fit with prism 9.
Y0 = 0.9981;
Plateau = 0.4970;
K = 0.0002185;
Y_iswap=(Y0 - Plateau)*exp(-K*iswap_in) + Plateau*ones(size(iswap_in));

Y0_pauli = 0.998;
Plateau_pauli = 0.4975;
K_pauli = 1.989e-5;
Y_pauli=(Y0_pauli - Plateau_pauli)*exp(-K_pauli*pauli_in) + Plateau_pauli*ones(size(pauli_in));



figure
box on
hold on

fill([iswap_in;flipud(iswap_in)],[iswap_f-iswap_err;flipud(iswap_f+iswap_err)],[0.70 0.80 0.940],'linestyle','none'); %
s1=line(iswap_in,iswap_f,'LineStyle',':','Color',[0 0.4470 0.7410],'marker','d','linewidth',2); %
s2=plot(iswap_in,Y_iswap,'-','Color',[0.4660 0.6740 0.1880],'linewidth',2);

fill([pauli_in;flipud(pauli_in)],[pauli_f-pauli_err;flipud(pauli_f+pauli_err)],[0.9 0.80 0.70],'linestyle','none');
s3=line(pauli_in,pauli_f,'LineStyle',':','Color',[0.9290 0.6940 0.1250],'marker','x','linewidth',2);
s4=plot(pauli_in,Y_pauli,'-','Color',[0.500 0.4250 0.0980],'linewidth',2);

legend('','Exp iswap','Fitted iswap','','Exp pauli','Fitted pauli','location','southeast','FontSize',18);
xlabel('Size','fontsize',18);
ylabel('Pr','fontsize',18);
set(gca,'fontsize',12);
