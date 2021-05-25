%% description
% This script sorts ellipsoids by how parallel their longest axes are
%
% Authors: Shreyas Kousik
% Created: 19 May 2021
% Updated: 25 May 2021
clear ; clc ;
%% user parameters
% rng seed
% rng(100) ; 

% number of ellipsoids to consider
n_ell = 3 ;

%% automated from here
% randomly generate a bunch of ellipsoids
E_list = cell(1,n_ell) ;
V_list = cell(1,n_ell) ;
l_list = cell(1,n_ell) ;

for idx = 1:n_ell
    Q_idx = make_random_covariance_matrix(2) ;
    
    [V_idx, l_idx] = eig(Q_idx) ;
    E_list{idx} = Q_idx ;
    V_list{idx} = V_idx ;
    l_list{idx} = diag(l_idx) ;
end

%% compute similarities based on longest semi-axis
V_1 = get_longest_semiaxis(V_list{1},l_list{1}) ;
sim_list = ones(1,n_ell) ;

for idx = 1:n_ell    
    V_idx = get_longest_semiaxis(V_list{idx},l_list{idx}) ;
    sim_list(idx) = abs(sin(V_idx'*V_1)) ;
end

%% plotting setup
th_vec = linspace(0,2*pi,200) ;
C = [cos(th_vec) ; sin(th_vec)] ;

%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ;

for idx = 1:n_ell
    % transform unit ball
    G = get_generators(E_list{idx}) ;
    P = G*C + [idx-1;0];
    
    % plot color based on d_list
    if idx > 1
        r = sim_list(idx) ;
        c = [(1-r) 0 r] ;
    else
        c = [0 0 1] ;
    end
    patch('faces',[1:200,1],'vertices',P','facecolor',c,'facealpha',1,...
        'edgealpha',0)
    
%     % plot arrow in largest eigenvector direction
%     V_idx = get_longest_semiaxis(V_list{idx},l_list{idx}) ;
%     plot_arrow([idx-1;0],[idx-1;0]+V_idx) ;
end

set_plot_linewidths(2) 

%% helper functions
function S = make_random_covariance_matrix(d)
% S = make_random_covariance_matrix(d)
%
% It does what it says on the box. The input d is the dimension, so the
% output is a d x d matrix.
%
% Authors: Shreyas Kousik
    Q = 2*rand(d) - 1;
    D = diag(rand(d,1)) ;
    S = Q*D*Q' ;
end

function G = get_generators(Q)
    G = inv(inv(Q)^(1/2)) ;
end

function v = get_longest_semiaxis(V,l)
    [~,l_max_idx] = max(l) ;
    v = V(:,l_max_idx) ;
end