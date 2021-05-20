%% description
% This script sorts ellipsoids by their volume and similarity of their
% "covariance" matrices
%
% Authors: Shreyas Kousik
% Created: 19 May 2021
% Updated: nah
clear ; clc ;
%% user parameters
% rng seed
rng(0) ; 

% number of ellipsoids to consider
n_ell = 4 ;

%% automated from here
% randomly generate a bunch of ellipsoids and compute their volumes
E_list = cell(1,n_ell) ;
v_list = nan(1,n_ell) ;
for idx = 1:n_ell
    E_list{idx} = rand_range(1,10).*make_random_covariance_matrix(2) ;
    v_list(idx) = vol_ellipsoid(E_list{idx}) ;
end

v_ratios = v_list./v_list(1) ;

%% compute ellipsoid similarities
d_list = nan(1,n_ell) ;

for idx = 1:n_ell
    d_list(idx) = KL_ellipsoid(E_list{1},E_list{idx}) ;
end

% d_list = d_list./v_ratios ;

d_list(1) = 0 ;
d_max = max(d_list) ;

%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ;

for idx = 1:n_ell
    E = ellipsoid(E_list{idx}) ;
    
    % plot color based on d_list
    if idx > 1
        r = d_list(idx)./d_max ;
        c = [0 r r] ;
    else
        c = [0 1 0] ;
    end
    h = plot(E) ;
    h.Color = c ;
end

set_plot_linewidths(1.5) 

%% helper functions
function S = make_random_covariance_matrix(d)
% S = make_random_covariance_matrix(d)
%
% It does what it says on the box. The input d is the dimension, so the
% output is a d x d matrix.
%
% Authors: Shreyas Kousik
    Q = rand(d);
    D = diag(rand(d,1)) ;
    S = Q*D*Q' ;
end

function v = vol_ellipsoid(Q)
    n = size(Q,1) ;
    v_n = (pi^(n/2))/gamma(n/2 + 1) ;
    v = v_n.*(det(inv(Q))^(1/2)) ;
end

function d = KL_ellipsoid(Q_1,Q_2)
    k = size(Q_1,1) ;
    d = (1/2)*(trace(inv(Q_1)*Q_2) - k + log(det(Q_1)/det(Q_2))) ;
end