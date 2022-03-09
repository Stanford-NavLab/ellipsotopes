%% description
% This script tests the lift-and-reduce strategy from the constrained
% zonotope paper [1]. In particular we leverage the fact that
% 2-ellipsotopes can be order-reduced in a nice way.
%
% [1] Scott, J.K., Raimondo, D.M., Marseglia, G.R. and Braatz, R.D., 2016.
%     Constrained zonotopes: A new tool for set-based estimation and fault
%     detection. Automatica, 69, pp.126-136.
%
% Authors: Shreyas Kousik
% Created: 3 Mar 2022
% Updated: 4 Mar 2022

%% user parameters
% rng seed
rng(0) ;

% original tope
n_dim = 2 ;
n_gen = 30 ;
n_con = 3 ;
n_I = 3 ;

%% automated from here
% make original etope
[E,c,G,A,b,I] = make_random_ellipsotope(2,n_dim,n_gen,n_con,n_I) ;

% lift
Gl = [G ; A] ;
cl = [c ; -b] ;
El = ellipsotope(2,cl,Gl,[],[],I) ;
nl = size(Gl,1) ; % dimension of lifted etope

% reduce by iterating through I
G_rdc = [] ;
I_rdc = {} ;
for idx = 1:n_I
    J = I{idx} ;
    G_idx = Gl(:,J) ;
    [nr,nc] = size(G_idx) ;
    if nc > nr
        G_idx = reduce_2_etope_generator_matrix(G_idx) ;
        I_idx = (1:nl) + get_max_index(I_rdc) ;
    else
        I_idx = (1:nc) + get_max_index(I_rdc) ;
    end
    G_rdc = [G_rdc, G_idx] ; 
    I_rdc = [I_rdc, {I_idx}] ;
end

% make new 'tope
c_rdc = cl(1:n_dim) ;
A_rdc = G_rdc((n_dim+1):end,:) ;
b_rdc = -cl((n_dim+1):end) ;
G_rdc = G_rdc(1:n_dim,:) ;
E_rdc = ellipsotope(2,c_rdc,G_rdc,A_rdc,b_rdc,I_rdc) ;

%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ;
plot(E)
plot(E_rdc,'color','r','linestyle','--','linewidth',3)

make_plot_pretty()