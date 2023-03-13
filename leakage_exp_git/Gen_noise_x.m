function rho_new = Gen_noise_x(rho, p, Nq)
%UNTITLED3 generate depolarize noise with pro p and leakage L
%the value for the specific elements: [1,2]:p_1,[2,1]:p_2,[2,2]:p_3
%p: with pro p, rho --> x rho x, where x=[0 0 0;0 0 1;0 1 0].
%
MAX=2^Nq;

max3q = 3^Nq;

leak_x = [1 0 0;0 0 1;0 1 0];

%generate all  x_i
leak_all = zeros(max3q,max3q,Nq);
for j = 1 : Nq
    leak_all(:,:,j) = sparse(eye(max3q));
end

for j = 1 : Nq
   leak_all(:,:,j) = kron(eye(3^(j-1)),kron(leak_x,eye(3^(Nq-j)))); 
   leak_all(:,:,j) = sparse(leak_all(:,:,j));
end

% com_project = eye(3^Nq);
% single_project = [1 0 0;0 1 0; 0 0 0];
% for j = 1 : Nq
%     com_project = kron(eye(3^(j - 1)), kron(single_project, eye(3^(Nq - j))));
% end

alp_0 = 2 - (1+p)^Nq;
rho_new = sparse(alp_0 * rho);
% rho_new = alp_0 * ((1 - p_dep) * sparse(rho) + p_dep/2^Nq * sparse(com_project));
% alp = 1;
% display(rho_new);

rho = sparse(rho);


all_bit = (1:MAX-1);
% all_gray = bin2gray(all_bit, 'qam', MAX);
all_string = dec2bin(all_bit)-'0';
for j = 1 : MAX-1
    rho_cur = sparse(eye(max3q));
   for k = 1 : Nq
       if all_string(j,k) == 1
            rho_cur = leak_all(:,:,k) * rho * leak_all(:,:,k)';
       end
   end
   alp = sum(all_string(j,:));
%    fprintf('alp=%d, coe =%f, ',alp, p^alp);
   rho_new = rho_new + rho_cur*p^alp;
%    fprintf('j=%d, rho_new = ',j);
%    display(rho_new);
%    display(rho_cur);

end



end

