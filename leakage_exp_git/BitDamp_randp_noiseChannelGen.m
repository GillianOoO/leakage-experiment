function [pp,p, ErrorChannel] = BitDamp_randp_noiseChannelGen(Nq, ExpN,noise_strength,precision)
%UNTITLED3 generate noise parameters 3^n by 3^n numbers in total.
% noise: noise_strength = [1e-3 (1 pm precision),1e-6(1 pm precision)];

pp = zeros(ExpN, ExpN);
p = zeros(2,Nq); %p(1:i) = pi and p(2:i)=qi

ErrorChannel = zeros(ExpN, ExpN, ExpN, ExpN);
%%%ErrorChannel(:, :, ExpN, ExpN) be E_0
%%%ErrorChannel(:, :, j, k) be the E_jk


for k = 1 : ExpN
    for j = 1 : ExpN
        if j == k
            continue;
        end
        for i = 1 : Nq
            if pp(k,j) ~= 0
                break;
            end
            if (In_compute_num(k-1) == 1 && If_IthqubitL(j-1,i,Nq)==1) ||...
                 (If_IthqubitL(k-1,i,Nq) == 1 && In_compute_num(j-1) == 1)    
                pp(k,j) = noise_strength(2)*(1 + rand()*precision);
                flag = 2;
                if In_compute_num(k-1) == 1 && If_IthqubitL(j-1,i,Nq)==1
                    flag = 1;
                end
                p(flag, i) = p(flag, i) + pp(k,j);
                ErrorChannel(j,k,k,j) = sqrt(pp(k,j));
            elseif (In_compute_num(k-1) && In_compute_num(j-1))|| ...
                    (If_IthqubitL(k-1,i,Nq)==1 && If_IthqubitL(j-1,i,Nq)==1)
                pp(k,j) = noise_strength(1)*(1 + rand()*precision);
                ErrorChannel(j,k,k,j) = sqrt(pp(k,j));
            end

        end

    end
end

%generate E0
for k = 1 : ExpN
    temp = 1 - sum(pp(k,:));
    ErrorChannel(k,k,ExpN, ExpN) = sqrt(temp);
end


end