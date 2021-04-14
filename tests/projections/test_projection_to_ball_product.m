%% description
% This script tests driving points towards the boundary of a ball defined
% by both the p-norm and the index set of a (generalized) ellipsotope.
%
% See also: test_projection_to_constraint_and_ball.m
%
% Authors: Shreyas Kousik
% Created: 4 Apr 2021
% Udpated: 5 Apr 2021
clear;clc
%% user parameters
% rng seed
rng(0) ;

% p norm
p_norm = 2 ;

% index set
I = {[1,2],[3]} ;

% number of points
n_P = 5000 ;

%% automated from here
% sanity check the index set
[I_chk,n_I,n_dim] = check_index_set_validity(I) ;
if ~I_chk
    error(['The index set is not valid! ',...
        'It should contain each dimension of the hyperball product ',...
        'exactly once.'])
end

% create a bunch of points in the space
P_orig = 2*rand(n_dim,n_P) - 1 ;

%% construct points in generalized e'tope ball
% make random list of which norm to enforce per point
idxs_J_to_enf = rand_int(1,n_I,[],[],1,n_P) ;

% project points to surface of product of balls
tic
P_proj = project_points_to_ball_product(P_orig,p_norm,I,idxs_J_to_enf) ;
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