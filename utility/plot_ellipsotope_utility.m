function plot_ellipsotope_utility(p_norm,c,G,A,b,I,varargin)
% plot_ellipsotope_utility(p,c,G,A,b,I)
%
% This function just reimplements the code in the test_plot_ellipsotope.m
% script for testing purposes.
%
% Authors: Shreyas Kousik
% Created: 13 Apr 2021
% Updated: nah


[I_chk,~,n_dim] = check_index_set_validity(I) ;
if ~I_chk
    error(['The index set is not valid! ',...
        'It should contain each dimension of the hyperball product ',...
        'exactly once.'])
end

% get number of constraints
if ~isempty(A)
    [~,n_dim_con] = size(A) ;
    
    % sanity check the constraints
    if n_dim_con ~= n_dim
        error(['The constraints defined by (A,b) must have as many columns ',...
            'as there are dimensions in the index set (i.e., the hyperball ',...
            'space of the ellipsotope.'])
    end
end

% create a bunch of points in the ball space
n_P = 1000 ;
P = 2*rand(n_dim,n_P) - 1 ;

% project points to intersection of ball product and linear subspace
[P,n_P] = project_points_to_ball_product_and_linear_subspace(P,p_norm,A,b,I) ;

% map points to workspace
P = repmat(c,1,n_P) + G*P ;
K = convhull(P') ;
P = P(:,K) ;
P = points_to_CCW(P) ;
n_P = size(P,2) ;

% plot ellipsotope
patch('faces',[1:n_P,1],'vertices',P','facecolor','b','facealpha',0.1,...
    'linewidth',1.5,'edgecolor','b',varargin{:})

end