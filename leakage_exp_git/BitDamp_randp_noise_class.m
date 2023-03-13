function [p,q]=BitDamp_randp_noise_class(pp, Nq)
%generate p and q with pp

ExpNq = length(pp);

p = zeros(1,Nq);
q = zeros(1,Nq);
for k = 1 : ExpNq
    for j = 1 : ExpNq
        if k == j
            continue;
        end
        for i = 1 : Nq
            if In_compute_num(k) == 1 && If_IthqubitL(j,i,Nq)==1
                p(i) = p(i) + pp(k,j);
            end
            if If_IthqubitL(k,i,Nq) == 1 && In_compute_num(j) == 1
                p(i) = p(i) + pp(j,k);
            end
        end
    end

end


end