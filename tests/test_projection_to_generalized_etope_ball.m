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
% p norm
p_norm = 2 ;

% index set
I = {[1,2],[3,4]} ;

% number of points
n_P = 5000 ;

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
P_proj = nan(size(P_orig)) ;

%% construct points that obey projection
tic
for idx_P = 1:n_P
    % for each point (we can vectorize this to make it faster)...
    p_idx = P_orig(:,idx_P) ;
    
    % pick at random which
    idx_J_to_enf = rand_int(1,n_I) ;
    
    for idx_I = 1:n_I
        J = I{idx_I} ;
        
        v = vecnorm(p_idx(J),p_norm) ;
        
        if idx_I == idx_J_to_enf
            % enforce the norm
            p_idx(J) = p_idx(J)./v ;
        else
            % check if the norm needs to be enforced
            if v > 1
                p_idx(J) = p_idx(J)./v ;
            end
        end
    end
    
    P_proj(:,idx_P) = p_idx ;
end
toc

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