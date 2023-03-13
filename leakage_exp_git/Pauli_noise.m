clear all;

syms p1 p2 p3
p3 = 0;
p1=0;

noise = [1 - p1-p3, p3, p1;
    p3, 1 - p2 - p3, p2;
    p1, p2, 1- p1 - p2];

P = [1/2 1/2 0;1/2 1/2 0;0 0 1];


F_mat = P * noise;

[states, vals] = eig(F_mat);