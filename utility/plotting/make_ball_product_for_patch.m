function [F,V,V_all] = make_ball_product_for_patch(p_norm,I,proj_dims,n_V)
% [F,V] = make_ball_product_for_patch(p_norm,I)
% [F,V] = make_ball_product_for_patch(p_norm,I,proj_dims)
% [F,V] = make_ball_product_for_patch(p_norm,I,proj_dims,n_V)
% [F,V,V_all] = make_ball_product_for_patch(...)
%
% Make faces F and vertices V to be used with patch for the ball product
% defined by the given p-norm and index set I. Optional inputs are the
% projection dimensions (1:3 by default) and the number of points to use in
% creating the patch object (10000 by default). The optional third output
% is the full-dimensional set of vertices on the boundary of the ball
% product.
%
% See also: project_points_to_ball_product
% 
% Authors: Shreyas Kousik
% Created: 6 Apr 2021
% Updated: nuu

    % set default inputs
    if nargin < 3
        proj_dims = 1:3 ;
    end
    
    if nargin < 4
        n_V = 10000 ;
    end
    
    % get dimensions of things
    n_dim = max(cell2mat(I)) ;
    n_I = length(I) ;
    
    % generate a bunch of points in the n_dim unit hypercube
    V = 2*rand(n_dim,n_V) - 1 ;
    
    % for each point pick a random index subset to project such that the
    % point is now on the boundary of the ball product
    idxs_J_V_E = rand_int(1,n_I,[],[],1,n_V) ;
    
    % project all the points to the boundary
    V_all = project_points_to_ball_product(V,p_norm,I,idxs_J_V_E) ;
    
    % get the projection dimensions and transpose the output
    V = V_all(proj_dims,:)' ;
    
    % get the convex hull of the resulting points
    F = convhull(V) ;
end