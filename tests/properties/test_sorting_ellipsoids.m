%% description
% This script sorts ellipsoids by the similarity metrics in [1]
%
% [1] Moshtaghi, M., Havens, T.C., Bezdek, J.C., Park, L., Leckie, C.,
%     Rajasegarar, S., Keller, J.M. and Palaniswami, M., 2011. Clustering
%     ellipses for anomaly detection. Pattern Recognition, 44(1), pp.55-69.
%
% Authors: Shreyas Kousik
% Created: 19 May 2021
% Updated: 21 May 2021
clear ; clc ;
%% user parameters
% rng seed
% rng(100) ; 

% number of ellipsoids to consider
n_ell = 3 ;

%% automated from here
% randomly generate a bunch of ellipsoids
E_list = cell(1,n_ell) ;
R_list = cell(1,n_ell) ;
l_list = cell(1,n_ell) ;

for idx = 1:n_ell
    Q_idx = make_random_covariance_matrix(2) ;
    
    [R_idx, l_idx] = eig(Q_idx) ;
    E_list{idx} = Q_idx ;
    R_list{idx} = R_idx ;
    l_list{idx} = diag(l_idx) ;
end

%% compute similarities to first ellipsoid
R_1 = R_list{1} ;
l_1 = l_list{1} ;

sim_list = ones(1,n_ell) ;

for idx = 1:n_ell
%     % compute orientation similarity
    R_idx = R_list{idx} ;
    th_1 = acos(diag(R_1'*R_idx)) ;
    th_2 = acos(diag(R_idx'*R_1)) ;
    val = min(vecnorm(sin([th_1,th_2]),2)) ;
    sim_idx = exp(-val) ;
    sim_list(idx) = sim_idx ;
end

%% compute similarities based on longest semi-axis
G_1 = get_generators(E_list{1}) ;

%% sort by similarity
[sim_list,idxs] = sort(sim_list,'descend') ;
E_list = E_list(idxs) ;

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
end

plot_arrow(zeros(2,1),R_1(:,1))
plot_arrow(zeros(2,1),R_1(:,2))
plot_arrow(zeros(2,1),G_1(:,1))
plot_arrow(zeros(2,1),G_1(:,2))

set_plot_linewidths(2) 

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

function G = get_generators(Q)
    G = inv(inv(Q)^(1/2)) ;
end