%% description
% This script demonstrates a (rough) metric for estimating the size of
% an ellipsotope (vs. a zonotope and an ellipsoid)
%
% Authors: Shreyas Kousik
% Created: 22 Dec 2021
% Updated: 23 Dec 2021
clear ; clc ;
%% user parameters
% random number generator seed
rng(0)

% 'tope specs
p_norm = 2 ;
n_dim = 2 ;
n_gen = 4 ;
n_con = 2 ;

% ray tracing
n_ray = 10 ;

%% automated from here
% make etope
[E,c,G,A,b,I] = make_random_ellipsotope(p_norm,n_dim,n_gen,n_con) ;

% make zonotope and ellipsoid
I_Z = num2cell(1:n_gen) ;
I_L = {1:n_gen} ;
Z = ellipsotope(p_norm,c,G,A,b,I_Z) ;
L = ellipsotope(p_norm,c,G,A,b,I_L) ;

%% estimate sizes
% estimate ellipsotope size
[v_E,P_E,~,D,p_out,x_out] = estimate_size_ray_tracing_ND(E,200) ;
disp(['Ellipsotope size: ',num2str(v_E,'%0.2f')])

% estimate zonotope size
[v_Z,P_Z] = estimate_size_ray_tracing_ND(Z,D,p_out,x_out) ;
disp(['Zonotope size: ',num2str(v_Z,'%0.2f')])

% estimate ellipsoid size
[v_L,P_L] = estimate_size_ray_tracing_ND(L,D,p_out,x_out) ;
disp(['Ellipsoid size: ',num2str(v_L,'%0.2f')])


%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ;

plot(E,'color','b')
plot(Z,'color','r')
plot(L,'color','g')

plot_path(P_E,'b.')
plot_path(P_Z,'r.')
plot_path(P_L,'g.')

make_plot_pretty()