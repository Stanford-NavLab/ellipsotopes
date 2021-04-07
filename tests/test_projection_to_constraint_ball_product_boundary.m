%% description
% This script tests solving for a scalar to push points to the boundary of
% a ball product.
%
% As of 9 PM 6 Apr 2021, this almost works!!
%
% Authors: Shreyas Kousik
% Created: 6 Apr 2021
clear;clc
%% user parameters
% rng seed
rng(0) ;

% p norm
p_norm = 2 ;

% number of points to create
n_P = 1 ;

% index set
I = {[1,2],3} ;

% constraint
% A = [1 1 1 ;
%     -1 -1 1] ;
% b = 0.1*ones(2,1) ;
A = [1 1 1] ;
b = 0.1 ;

%% automated from here
% sanity check the index set
[I_chk,n_I,n_dim] = check_index_set_validity(I) ;
if ~I_chk
    error(['The index set is not valid! ',...
        'It should contain each dimension of the hyperball product ',...
        'exactly once.'])
end

% get number of constraints
n_con = size(A,1) ;

%% project points to intersection of constraint and ball product
% get nullspace of constraint
K_con = null(A) ;
t_con = A\b ;
n_null_dim = size(K_con,2) ;

% create some random in the nullspace of the constraint
P_dir = K_con*(2*rand(n_null_dim,n_P) - 1) ;

% for each random vector, solve for the alpha that puts it to the boundary

%% testing solving for coefficients
u = P_dir ;
pows = 0:p_norm ;
coef_u = u.^pows(end:-1:1) ;
coef_t = t_con.^pows ;
coef_binoms = repmat(vec_nchoosek(p_norm),n_dim,1) ;
coef_poly = coef_binoms.*coef_u.*coef_t ;
coef_poly(:,end) = coef_poly(:,end) - 1 ;

coef_sols = arrayfun(@(idx) roots(coef_poly(idx,:)), 1:n_dim,'uniformoutput',false) ;
coef_sols = cell2mat(coef_sols) ;
coef_val = min(coef_sols(1,:)) ;

% create point
p_test = u*coef_sols(:)' + t_con ;

% test that at least one point obeys the norm thing
vecnorm_ball_product(p_test,p_norm,I)

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
plot_arrow(t_con,t_con + u)

% plot point
plot_path(p_test,'rx','markersize',10)

% plot hyperplanes
for idx = 1:n_con
    % get nullspace info
    K = K_each{idx} ;
    t = t_each{idx} ;
    
    % create hyperplane
    F_H = [1 2 3 4 1] ;
    V_H = 2.*[K' ; -K'] + repmat(t',4,1) ;
    patch('faces',F_H,'vertices',V_H,'facealpha',0.1','edgealpha',0,'facecolor','r') ;
end