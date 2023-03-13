 clear all;

 Nq = 3;
 ExpN = 3^Nq;
 rho = zeros(ExpN, ExpN);
 rho(1,1) = 1;
noise_strength = [1e-3,1e-6];% [1e-4 (1 pm precision),1e-7(1 pm precision)];
precision = 5e-1;
%get pp in 3^Nq * 3^Nq, p_i in 2 * Nq, ErrorChannel(:,:,i,j) denotes Eij,
%and ErrorChannel(:,:,ExpN,ExpNq)
[pp,p, ErrorChannel] = BitDamp_randp_noiseChannelGen(Nq, ExpN,noise_strength,precision);

save('pp_test_new.txt','pp', '-ascii');
save('p_test_new.txt','p', '-ascii');

[rho] = BitDampNoise_randp(rho, ErrorChannel, pp);
