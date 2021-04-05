%% description
% This script tests driving points towards the boundary of a ball defined
% by both the p-norm and the index set of a (generalized) ellipsotope.
%
% See also: test_projection_to_constraint_and_ball.m
%
% Authors: Shreyas Kousik
% Created: 4 Apr 2021
% Udpated: not yet
clear;clc
%% user parameters
% rng seed
rng(0) ;

% p norm
p_norm = 2 ;

% index set
I = {[1,2],[3,4],[5,6],[7,8,9]} ;

% number of points
n_P = 10000 ;

%% automated from here
% sanity check the index set
n_I = length(I) ;
D = [] ;
for idx = 1:length(I)
    D = [D, I{idx}] ;
end

% get the dimension
n_dim = max(D) ;

if length(unique(D)) < length(D) || length(D) < n_dim
    error('The index set is not valid!')
end

% create a bunch of points in the space
P_orig = 2*rand(n_dim,n_P) - 1 ;
P_proj = P_orig ;

%% construct points in generalized e'tope ball
tic

% make random list of which norm to enforce per point
idxs_J_to_enf = rand_int(1,n_I,[],[],1,n_P) ;

% iterate through the index list and enforce things on points...
for idx_J = 1:n_I
    % get current index from the index set
    J = I{idx_J} ;
    
    % get the norm of all points for the current index set
    V = vecnorm(P_proj(J,:),p_norm,1) ;
    
    % get all points on which to enforce the norm for this index, which
    % includes all points for which the norm is violated by being too big
    idxs_J = (idxs_J_to_enf == idx_J) | (V > 1) ;
    
    % enforce the norm
    P_proj(J,idxs_J) = P_proj(J,idxs_J) ./ repmat(V(idxs_J),length(J),1) ;   
end
toc

%% validate that all points obey the norms/indices

%% plotting
figure(1) ; clf ;

if n_dim == 3
    axis equal ; hold on ; view(3) ;
    
    % plot original points
    plot_path(P_orig,'r.','markersize',4)
    
    % plot projected points
    plot_path(P_proj,'b.','markersize',4)
    
    % labeling
    xlabel('x_1') ;
    ylabel('x_2') ;
    zlabel('x_3') ;
    legend('original points','points on boundary')
    
    set(gca,'fontsize',15)
else
    % create a bunch of 2-D subplots
    n_plot = n_dim - 1 ;
    for idx_plot = 1:n_plot
        % get the subplot
        subplot(n_plot,1,idx_plot) ; axis equal ; hold on ; grid on ;
        
        % get the current projection dimensions
        proj_dims = [idx_plot,idx_plot+1] ;
        
        % plot original points
        plot_path(P_orig(proj_dims,:),'r.','markersize',4)
        
        % plot projected points
        plot_path(P_proj(proj_dims,:),'b.','markersize',4)
        
        % labeling
        xlabel(['x_',num2str(proj_dims(1))])
        ylabel(['x_',num2str(proj_dims(2))])
        legend('original points','points on boundary')
        
        set(gca,'fontsize',15)
    end
end