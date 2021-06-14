%% description
clear ; clc
%% user parameters
% ellipsotope definitions
c_1 = [0;0] ;
G_1 = [+1 -1 ;
       -1 +0] ;
   
c_2 = [0;0] ;
G_2 = [-1 -1 ;
       -1 +2] ;
   
c_3 = [0;0] ;
G_3 = [-2 -1 ;
       -1 +1] ;

%% automated from here
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
E_prime = E_left_out + E_rdc ;

%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ;

plot(E_1,'color','r')
plot(E_2,'color','b')
plot(E_3,'color','g')
plot(E,'color',0.5*ones(1,3))
plot(E_prime,'color',[0.5 0.5 0])