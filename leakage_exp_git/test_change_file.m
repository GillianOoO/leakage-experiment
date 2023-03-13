clear all;

f_in = 'Fidelity_ILRB_pauli_original.txt';
f_out = 'Fidelity_ILRB_pauli.txt';

data = load(f_in);

data_out = data(1:4:400,:);

save(f_out,'data_out','-ascii');