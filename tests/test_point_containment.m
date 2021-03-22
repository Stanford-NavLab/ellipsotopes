% Test point containment checking for ellipsotopes

c = [0;0];
G = [2 0 0 0;
     0 1 0 0];
A = [2 0 0 0;
     0 1 0 0];
b = [1; 0];
I = {[1,2],[3,4]};

E = ellipsotope(2,c,G,A,b,I);

% sample some points
n_P = 50;
b_samp = [-2,2];
P = make_grid(repmat(b_samp,1,2),n_P*ones(1,2));

% test the points
H = [G; A];
N = null(H);
for i = length(P)
    x = P(:,i);
    h = [x - c; b];
    