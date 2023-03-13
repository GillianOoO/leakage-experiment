clear all;
%%generate 2n random number in [A,B]

Nq= 5;
A = 10^-4;
B = 1.5 * 10^-4;

ErrVal = zeros(2,Nq);
for j = 1 : 2
    for k = 1 : Nq
        ErrVal(j,k) = A + rand()*(B-A); %in A~A+B
    end
end

save('ErrVal.txt','ErrVal','-ascii');

