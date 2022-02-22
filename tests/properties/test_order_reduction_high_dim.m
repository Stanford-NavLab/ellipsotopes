%% description
% This script creates a figure illustrating the order reduction heuristic
% for 2-ellipsotopes in higher dimensions.
%
% See also: figure_order_reduction_2_etopes.m
%
% Authors: Shreyas Kousik
% Created: 21 Feb 2022
% Updated: --
clear ; clc
%% user parameters
% random number generator seed
% rng(0) ;

% number of topes to test
n_topes = 12 ;

% dimension
n_dim = 10 ;

% number of rays to use for volume estimation
n_ray = 200 ;

%% automated from here
% create initial ellipsotope
c = zeros(n_dim,1) ;
G = inv(sqrt(make_random_covariance_matrix(n_dim))) ;
G = G./max(G(:)) ;

% create all ellipsotopes
E_cell = cell(1,n_topes) ;
E_cell{1} = ellipsotope(2,c,G) ;

for idx = 2:n_topes
    R_idx = RandOrthMat(n_dim) ; % make_random_orthonormal_matrix(n_dim) ;
    s_idx = 1 ; %(2*rand(1) - 1) + 1 ;
    G_idx = s_idx.*R_idx*G ;
    E_cell{idx} = ellipsotope(2,c,G_idx) ;
end

%% creating a new heuristic
% possible combinations of topes
combs = combinator(n_topes,2,'c') ;
n_combs = size(combs,1) ;

% set up ray directions for estimating volume of 'topes
D = rand(n_dim,n_ray) ;

% for each combinations...
vols_MVOE = nan(1,n_combs) ; % volume of MVOE
vols_rdc = nan(1,n_combs) ; % volume of reduced 'tope
vols_heur = nan(1,n_combs) ; % heuristic value
bts = nan(1,n_combs) ;
tic
for idx = 1:n_combs
    % get ellipsoids to compare
    idx_i = combs(idx,1) ;
    idx_j = combs(idx,2) ;
    
    % get the centers and generator matrices
    c_i = E_cell{idx_i}.center ;
    c_j = E_cell{idx_i}.center ;
    G_i = E_cell{idx_i}.generators ;
    G_j = E_cell{idx_j}.generators ;
    
    % compute the MVOE
    [G_rdc,~,bt] = make_MVOE_generator_matrix(G_i,G_j) ;
    vols_MVOE(idx) = ellipsoid_volume_from_generator_matrix(G_rdc) ;
    bts(idx) = bt ;
    
    % compute heuristic
    % get shape matrices AS IN THE HALDER PAPER! (inverse of etope paper)
    Q_i = inv(pinv(G_i)'*pinv(G_i)) ;
    Q_j = inv(pinv(G_j)'*pinv(G_j)) ;
    Q_MVOE = 2*(Q_i + Q_j) ; % assume bt = 1 ;
    vols_heur(idx) = det(inv(Q_MVOE)) ;
    
    
    % create Minkowski sum 'tope
    idxs = 1:n_topes ; % indices of all topes
    idxs(combs(idx,:)) = [] ;
    E_temp = ellipsotope(2,c_i + c_j,G_rdc) ;
    for idx_inner = idxs
        E_temp = E_temp + E_cell{idx_inner} ;
    end
    
    % estimate volume as length along ray directions
    n_gen = E_temp.n_generators ;
%     vols_rdc(idx) = E_temp.estimate_size_ray_tracing_ND(D,zeros(n_dim,1),zeros(n_gen,1)) ;
%     toc
end

%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ;

plot(vols_MVOE,1./vols_heur,'bx')
make_plot_pretty()