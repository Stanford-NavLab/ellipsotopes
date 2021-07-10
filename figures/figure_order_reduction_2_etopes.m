%% description
% This script creates a figure illustrating the order reduction heuristic
% for 2-ellipsotopes.
%
% Authors: Shreyas Kousik
% Created: 13 June 2021
% Updated: 9 July 2021
clear ; clc
%% user parameters
% ellipsotope definitions
c_1 = [0;0] ;
G_1 = [+3 +2 ;
        -1 +0] ;
   
c_2 = [0;0] ;
% G_2 = [-1 -1 ;
%        -1 +6] ;
%    
c_3 = [0;0] ;
% G_3 = [-7 -1 ;
%        -1 +1] ;

%% automated from here
% rotate the generators of E_1
G_2 = rotation_matrix_2D(pi/6)*G_1 ;
G_3 = rotation_matrix_2D(pi/2)*G_1 ;

% make component etopes
E_1 = ellipsotope(2,c_1,G_1) ;
E_2 = ellipsotope(2,c_2,G_2) ;
E_3 = ellipsotope(2,c_3,G_3) ;

% minkowski sum all of them
E = E_1 + E_2 + E_3 ;

%% compute order reduction statistic
% gather the cs and Gs and Es
c = {c_1,c_2,c_3} ;
G = {G_1,G_2,G_3} ;
E_cell = {E_1,E_2,E_3} ;

% compute ellipsoid axes for each tope
V = cell(1,3) ;
for idx = 1:3
    G_idx = G{idx} ;
    G_inv = pinv(G_idx) ;
    Q = G_inv'*G_inv ;
    [V_idx,E_idx] = eig(inv(Q)) ;
    
    % get axis for larest eigenvalue
    E_idx = diag(E_idx) ;
    [~,idx_max] = max(E_idx) ;
    
    V{idx} = V_idx(:,idx_max) ;
end

% create all (i,j) pairs
combs = combinator(3,2,'c') ;
n_idx = size(combs,1) ;

% iterate over pairs to compute statistic
rdc_heur = nan(3,1) ;
for idx = 1:n_idx
    % get ellipsoids to compare
    ij = combs(idx,:) ;
    
    % get longest axes
    v_i = V{ij(1)} ;
    v_j = V{ij(2)} ;
    
    % compute heuristic
    rdc_heur(idx) = abs(v_i'*v_j) ;
end

[~,idx_max] = max(rdc_heur) ;
ij = combs(idx_max,:) ;
idx_i = ij(1) ;
idx_j = ij(2) ;

% find which tope was left out
idxs = 1:3 ;
log_i = idxs == idx_i ;
log_j = idxs == idx_j ;
idx_left_out = idxs(~(log_i | log_j)) ;

%% compute reduced etope
% compute MVOE generator matrix for max of heuristic
G_rdc = make_MVOE_generator_matrix(G{idx_i},G{idx_j}) ;

% create center
c_rdc = c{idx_i} + c{idx_j} ;

% get MVOE as an etope
E_rdc = ellipsotope(2,c_rdc,G_rdc) ;

% compute reduced etope
E_left_out = E_cell{idx_left_out} ;
E_tilde = E_left_out + E_rdc ;

%% compute other reduction possibilities
E_possible = cell(1,3) ;
vols = nan(1,3) ;
for idx = 1:3
    idxs = 1:3 ;
    idx_comb = combs(idx,:) ;
    log_i = idxs == idx_comb(1) ;
    log_j = idxs == idx_comb(2) ;
    idx_left_out_new = idxs(~(log_i | log_j)) ;
    
    G_rdc = make_MVOE_generator_matrix(G{idx_comb(1)},G{idx_comb(2)}) ;
    
    E_rdc_idx = ellipsotope(2,c_rdc,G_rdc) ;
    vols(idx) = det(G_rdc)*pi 
    E_left_out_idx = E_cell{idx_left_out_new} ;
    
    E_possible{idx} = E_left_out_idx + E_rdc_idx ;
end

%% plotting input topes
fh_a = figure(1) ; clf ; axis equal ; hold on ; grid on ;
h_3 = plot(E_3,'color',[0.7 0 0.7]) ;
h_2 = plot(E_2,'color',[0 0.7 0.7]) ;
h_1 = plot(E_1,'color',[0.5 0.5 0]) ;
legend([h_1 h_2 h_3],{'E_1','E_2','E_3'},'location','northwest')
axis tight
set_plot_fontsize(15) ;

%% plotting Mink sum figure
fh_b = figure(2) ; clf ; axis equal ; hold on ; grid on ;
for idx = 1:3
    plot(E_possible{idx},'color','r')
end
plot(E_tilde,'facecolor',[0.7 0.7 1],'edgecolor',[0 0 0.3],'facealpha',1,...
    'linewidth',3)
plot(E,'facealpha',0.0,'linestyle','--','edgecolor',[0 0 0.3],'linewidth',3,...
    'edgealpha',1) ;

axis tight
set_plot_fontsize(15) ;

%% save figures
save_figure_to_png(fh_a,'order_reduc_heur_2-etope_a.png')
save_figure_to_png(fh_b,'order_reduc_heur_2-etope_b.png')

%% compute areas of the shapes (only possible after plotting)
a = compute_approx_etope_area(E) ;
a_tilde = compute_approx_etope_area(E_tilde)/a
a_alt_1 = compute_approx_etope_area(E_possible{2})/a
a_alt_2 = compute_approx_etope_area(E_possible{3})/a

%% compute Hausdorff distances between the shapes
d_12 = compute_haus_dist(E,E_tilde)
d_13 = compute_haus_dist(E,E_possible{2})
d_23 = compute_haus_dist(E,E_possible{3})

%% helper functions
function a = compute_approx_etope_area(E)
    V = E.plot_handle.Vertices ;
    a = area(alphaShape(V)) ;
end

function d = compute_haus_dist(E_1,E_2)
% based on HausdorffDist from MathWorks FileExchange
    V_1 = E_1.plot_handle.Vertices ;
    V_2 = E_2.plot_handle.Vertices ;
    
    % get distances between points
    D = dist_points_to_points(V_1',V_2') ;
    
    % Obtain the value of the point, p, in P with the largest minimum distance
    % to any point in Q.
    vp = max(min(D,[],2));
    % Obtain the value of the point, q, in Q with the largets minimum distance
    % to any point in P.
    vq = max(min(D,[],1));

    d = max(vp,vq);
end
