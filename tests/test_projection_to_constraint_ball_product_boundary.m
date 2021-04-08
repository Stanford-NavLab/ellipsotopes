%% description
% This script tests solving for a scalar to push points to the boundary of
% the intersection of a linear subspace with a ball product.
%
% To do:
%   - expand the code to higher than 3-D ball products
%   - parallelize with respect to all the input points
%   - check for when the linear subspace is of dimension 1
%
% Authors: Shreyas Kousik
% Created: 6 Apr 2021
% Updated: 7 Apr 2021

clear;clc
%% user parameters
% rng seed
rng(0) ;

% p norm
p_norm = 2 ;

% number of points to create
n_P = 500 ;

% index set
I = {[1,2],3} ;

% constraint (it's prettier with just one constraint, but the code below
% works for an arbitrary number of constraints)
% A = [1 1 1 ;
%     -1 -1 1] ;
% b = 0.1*ones(2,1) ;
A = [-1 1 -1] ;
b = 0.2 ;

% tolerance for point being on boundary
tol_point_on_bdry = 1e-12 ;

%% automated from here
% sanity check the index set
[I_chk,n_I,n_dim] = check_index_set_validity(I) ;
if ~I_chk
    error(['The index set is not valid! ',...
        'It should contain each dimension of the hyperball product ',...
        'exactly once.'])
end

% get sizes of things
n_con = size(A,1) ;

%% project points to intersection of constraint and ball product
% get nullspace of constraint
K_con = null(A) ;
t_con = A\b ;
n_null_dim = size(K_con,2) ;

% create some random in the nullspace of the constraint
P_dir = K_con*(2*rand(n_null_dim,n_P) - 1) ;

%% testing solving for coefficients
tic
% for each random vector, solve for the coefficient that gets it to the
% boundary of the ball product and constraint intersection
P_out = [] ;

for idx_P = 1:n_P
    % get current point
    p_idx = P_dir(:,idx_P) ;
    
    % set 
    coef_sols = [] ;
    for idx_J = 1:n_I
        % get current index set
        J = I{idx_J} ;
        
        % get nullspace vectors for the current dimensions
        u = p_idx(J,:) ;
        
        % construct polynomial (binomial) coefficients given p-norm, since
        % we're solving sum((coef*u(J) - t(J))^p) = 1 for coef
        pows = 0:p_norm ;
        
        % apply binomial theorem
        coef_u = u.^pows(end:-1:1) ;
        coef_t = t_con(J).^pows ;
        coef_binoms = vec_nchoosek(p_norm) ; % binomial coefficients
        
        % sum the terms to get all the coefficients for the (scalar) coef poly
        coef_poly = coef_binoms.*sum(coef_u.*coef_t,1) ;
        
        % subtract 1 from both sides to get the polynomial in standard form
        coef_poly(:,end) = coef_poly(:,end) - 1 ;
        
        % solve for roots
        coef_sol = roots(coef_poly) ;
        
        % remove imaginary roots
        coef_sol_imag = imag(coef_sol) ;
        coef_sol(coef_sol_imag ~= 0) = [] ;
        
        % stack 'em
        coef_sols = [coef_sols, coef_sol(:)'] ;
    end
    
    % create candidate points
    p_test = p_idx*coef_sols + t_con ;
    
    % test to find which point obeys the norm
    v_idx = vecnorm_ball_product(p_test,p_norm,I) ;
    v_test = abs(v_idx - 1) < tol_point_on_bdry ;
    % idx_ok = find(v_test,1) ;
    
    % save the points that are ok
    P_out = [P_out, p_test(:,v_test)] ;
end
toc

%% plotting setup
% create ball product for patch
[F_E,V_E] = make_ball_product_for_patch(p_norm,I) ;

% get nullspace and translation of each linear constraint
K_each = cell(1,n_con) ;
t_each = cell(1,n_con) ;
for idx = 1:n_con
    K_each{idx} = null(A(idx,:)) ;
    t_each{idx} = A(idx,:)\b(idx) ;
end

%% plotting
figure(1) ; clf ; axis equal ; grid on ; hold on ; view(3) ;

% plot ball product
patch('faces',F_E,'vertices',V_E,'facealpha',0.1','edgealpha',0,'facecolor','b') ;

% plot t_con
plot_path(t_con,'b^','markerfacecolor','b')
% plot_arrow(t_con,t_con + P_dir)

% plot points on boundary
% plot_path(p_test,'rx','markersize',10)
plot_path(P_out,'b.','markersize',10)

% plot hyperplanes
for idx = 1:n_con
    % get nullspace info
    K = K_each{idx} ;
    t = t_each{idx} ;
    
    % create hyperplane
    F_H = [1 2 3 4 1] ;
    V_H = 3.*[K' ; -K'] + repmat(t',4,1) ;
    patch('faces',F_H,'vertices',V_H,'facealpha',0.1','edgealpha',0,'facecolor','r') ;
end