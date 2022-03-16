function [rdc_heur,rdc_pair_idxs,vols] = compute_order_reduction_heuristic(E_cell,flag_compute_volumes)
% rdc_heur = compute_order_reduction_heuristic(E_list)
%
% Given a cell array of basic ellipsotopes (i.e., ellipsoids) in E_list,
% compute the order reduction heuristic.%
%
% IT TURNS OUT THIS DOESN'T WORK SO WELL IN HIGHER DIMENSIONS -- OOPS!
%
% Authors: Shreyas Kousik
% Created: 21 Feb 2022
% Updated: 15 Mar 2022 (added warning)

warning('Hey! This heuristic is broken! Don''t use it!')

if nargin < 2
    flag_compute_volumes = false ;
end

% total number of topes
n_topes = length(E_cell) ;

% compute ellipsoid axes for each 'tope
V = cell(1,n_topes) ;
for idx = 1:n_topes
    G_idx = E_cell{idx}.generators ;
    G_inv = pinv(G_idx) ;
    Q = G_inv'*G_inv ;
    [V_idx,E_idx] = eig(inv(Q)) ;
    
    % get axis for larest eigenvalue
    E_idx = diag(E_idx) ;
    [~,idx_max] = max(E_idx) ;
    
    V{idx} = V_idx(:,idx_max) ;
end

% possible combinations of topes
combs = combinator(n_topes,2,'c') ;
n_combs = size(combs,1) ;

% set up to save volumes if necessary
if flag_compute_volumes
    vols = nan(1,n_combs) ;
end

% iterate over pairs to compute statistic
rdc_heur = nan(1,n_topes) ;
for idx = 1:n_combs
    % get ellipsoids to compare
    ij = combs(idx,:) ;
    
    % get longest axes
    v_i = V{ij(1)} ;
    v_j = V{ij(2)} ;
    
    % compute heuristic
    rdc_heur(idx) = abs(v_i'*v_j) ;
    
    if flag_compute_volumes
        G_rdc = make_MVOE_generator_matrix(E_cell{ij(1)}.generators,E_cell{ij(2)}.generators) ;
        vols(idx) = det(G_rdc)*pi ;
    end
end

[~,idx_max] = max(rdc_heur) ;
rdc_pair_idxs = combs(idx_max,:) ;
end