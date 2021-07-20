%% automated from here
rng(0)
[E,c,G,A,b,I] = make_random_ellipsotope(2,2,10,4,3) ;

%% block diagonalize A
% create Gm and Lm matrices
Gm = zeros(E.n_dimension,E.n_constraints) ;
Gm(1:2,1:2) = eye(2) ;
Lm = zeros(E.n_constraints,E.n_constraints) 
Lm(1) = 1

A - Lm*A

%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ;
plot(E)