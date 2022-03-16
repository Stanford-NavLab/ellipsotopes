%% description
% This script tests overapproximating an ellipsotope per [1, Prop. 5],
% which allows us to delete a constraint and a generator
%
% [1] Scott, J.K., Raimondo, D.M., Marseglia, G.R. and Braatz, R.D., 2016.
%     Constrained zonotopes: A new tool for set-based estimation and fault
%     detection. Automatica, 69, pp.126-136.
%
% Authors: Shreyas Kousik
% Created: 29 May 2021
% Updated: 16 Mar 2022 (using new function)
clear ; clc
%% user parameters
% rng
rng(0)

% etope specs
p_norm = 2 ;
n_dim = 2 ;
n_gen = 10 ;
n_con = 3 ;

%% automated from here
% make random etope
[E,c,G,A,b,I] = make_random_ellipsotope(p_norm,n_dim,n_gen,n_con) ;

% reduce
E_rdc = reduce_etope_constraint_and_generator(E) ;

%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ;

% plot reduced etope
plot(E_rdc,'facecolor','r','edgecolor','r','linestyle','--','facealpha',0.1)

% plot original etope
plot(E) ;