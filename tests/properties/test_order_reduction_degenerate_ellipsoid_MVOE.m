%% description
% This script demonstrates augmenting the generator matrices of an
% ellipsotope that has degenerate component ellipsoids (as may happen
% during lift-and-reduce order reduction). In particular, we make the
% degenerate ellipsoids no longer degenerate by augmenting them with
% additional random generators of very small magnitude, then compute the
% MVOE of the resulting component ellipsoids. It's just an idea for now.
%
% Authors: Shreyas Kousik
% Created: 14 Mar 2022
% Updated: -
clear ; clc ;
%% automated from here
G_i = 2*rand(4,2) - 1 ;
G_j = 2*rand(4,2) - 1 ;

% make G_i and G_j full rank
G_i_full = make_full_rank(G_i) ;
G_j_full = make_full_rank(G_j) ;


G_ij = make_MVOE_generator_matrix(G_i_full,G_j_full) ;

%% plotting setup
E = ellipsotope(2,zeros(4,1),[G_i,G_j],[],[],{[1 2],[3 4]}) ;
E_i = ellipsotope(2,zeros(4,1),G_i_full) ;
E_j = ellipsotope(2,zeros(4,1),G_j_full) ;
E_rdc = ellipsotope(2,zeros(4,1),G_ij) ;

%% plotting
figure(1) ; clf ;

plot(E,'color',[0 0 0],'plot_all_dims',true)
plot(E_rdc,'color',[1 0 0],'plot_all_dims',true)
plot(E_i,'color',[1 1 0],'plot_all_dims',true)
plot(E_j,'color',[0 1 1],'plot_all_dims',true)


%% helper functions
function G = make_full_rank(G,tol)
    [nr,nc] = size(G) ;
    if nargin < 2
        tol = 1e-6 ;
    end
    if nr > nc
        G = [G, tol.*rand(nr,nr-nc)] ;
    end
end