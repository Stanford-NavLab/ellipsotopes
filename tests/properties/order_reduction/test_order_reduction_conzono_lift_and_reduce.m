%% description
% Tests out the constrained zonotope paper's fancy nice order reduction
% methods (see [1]).
%
% [1] Scott, J.K., Raimondo, D.M., Marseglia, G.R. and Braatz, R.D., 2016.
%     Constrained zonotopes: A new tool for set-based estimation and fault
%     detection. Automatica, 69, pp.126-136.
%
% Authors: Shreyas Kousik
% Created: 15 Mar 2022
% Updated: --
clear ; clc
%% user parameters
% rng seed
rng(0)

% constrained zonotope specs
n_dim = 2 ;
n_gen = 15 ;
n_con = 5 ;

% number of generators to reduce
n_rdc = 3 ;

%% automated from here
% create the conzono as an ellipsotope
I = make_zonotope_index_set(n_gen) ;
Z = make_random_ellipsotope(2,n_dim,n_gen,n_con,I) ;

%% order reduction (see the appendix in [1])
% get propzzz
[p,c,G,A,b,I] = Z.get_properties() ;

% lift
G_l = [G ; A] ;
c_l = [c ; -b] ;
n_gen_l = size(G_l,2) ;
n_dim_l = length(c_l) ;

% reduce
G_rdc = reduce_zonotope_Chischi(G,n_rdc) ;

% reconstruct constrained zonotope
Z_rdc = ellipsotope(p,c_l,G_rdc,[],[],make_zonotope_index_set(size(G_rdc,2))) ;

% drop
Z_rdc = drop(Z_rdc,2) ;

%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ;

tic
plot(Z)
plot(Z_rdc,'color','r')
toc