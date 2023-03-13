clear all;

syms p1 p4 %pd2 pd3 %p1 p2 pd1 %p4 pd3


% p1 = 0;
p2 = 0;
pd1 = 0;
p3 = 2 * p2 + pd1;
%p4 = 2 * p1 + pd2;
pd3 = 0;

mat = [(4 - p3 - p4)/4, p3/4, p4/4, 0;
    p3/2, (2- p3 - pd3 - p1)/2, pd3/2, p1/2;
    p4/2, pd3/2, (2 - p4 - pd3)/2, p2/2;
    0, p1, p2, 1 - p1 - p2
    ];

[states, vals] = eig(mat);
vals = diag(vals);