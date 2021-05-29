%% description
% This script tests overapproximating an ellipsotope per [1, Prop. 5],
% which allows us to delete a constraint (wow!)
%
% [1] Scott, J.K., Raimondo, D.M., Marseglia, G.R. and Braatz, R.D., 2016.
%     Constrained zonotopes: A new tool for set-based estimation and fault
%     detection. Automatica, 69, pp.126-136.
%
% Authors: Shreyas Kousik
% Created: 29 May 2021
clear ; clc
%% user parameters
% rng
rng(0)

% etope specs
p_norm = 2 ;
n_dim = 2 ;
n_gen = 10 ;
n_con = 3 ;

%% autoamted from here
% make random etope
[E,c,G,A,b,I] = make_random_ellipsotope(p_norm,n_dim,n_gen,n_con) ;

% create Gamma and Lambda matrices
Gm = zeros(n_dim,n_con) ;
% Lm = 10*rand(n_con,n_con) - 5 ;
Lm = diag([1,zeros(1,n_con-1)]) ;

%% iterate over the constraints and eliminate one at a time
E_list = cell(1,n_con) ;

for idx = 1:n_con
% create Lambda
lm = zeros(n_con,1) ;
lm(idx) = 1 ;
Lm = diag(lm) ;
    
% create new etope
c_rdc = c + Gm*b ;
G_rdc = G - Gm*A ;
A_rdc = A - Lm*A ;
b_rdc = b - Lm*b ;

% delete first constraints
A_rdc = A_rdc(~lm,:) ;
b_rdc = b_rdc(~lm) ;

E_rdc = ellipsotope(p_norm,c_rdc,G_rdc,A_rdc,b_rdc,I) ;

E_list{idx} = E_rdc ;
end

%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ;

% plot all etopes
for idx = 1:n_con
    plot(E_list{idx},'facecolor','r','edgecolor','r','linestyle','--','facealpha',0.1)
end

% plot origetope
plot(E) ;