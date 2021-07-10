function [F,V] = make_patch_data_coeff_sampling(p,c,G,A,b,I,n_P)
% [F,V] = make_patch_data_coeff_sampling(p,c,G,A,b,I,n_P)
%
% This creates faces and vertices for a 2-D or 3-D ellipsotope by sampling
% in the coefficient space. The optional input argument n_P is the number
% of samples to generate, which is created by default based on the number
% of generators of the ellipsotope.
%
% Authors: Shreyas Kousik and Adam Dai
% Created: 20 May 2021
% Updated: 21 May 2021 (fixed bug in number of constraints if statement)

    % get sizes of things
    [n_dim,n_gen] = size(G) ;
    n_con = size(A,1) ;

    % sanity check
    if (n_dim ~= 2) && (n_dim ~= 3)
        error('This function only works for 2-D and 3-D ellipsotopes!')
    end
    
    % if n_P is not set, set it to default value
    if nargin < 7
        n_P = min(7^size(G,2),10^5);
    end
    
    % generate a bunch of points
    if n_con == 0
        switch n_gen
            case 2
                P = make_superellipse_2D(n_P,p);
            case 3
                P = make_superellipse_3D(n_P,p);
            otherwise
                P = make_unit_superellipse_ND(n_P,p,n_gen);
        end
    else
        % generate a bunch of random points 
        P = 2*rand(n_gen,n_P) - 1 ;
        
        % project points to ball and linear subspace boundary (the existing
        % function will handle index sets and empty linear subspaces properly,
        % it turns out)
        [P,n_P] = project_points_to_ball_product_and_linear_subspace(P,p,A,b,I) ;
    end
    
    % map points etope workspace
    P = repmat(c,1,n_P) + G*P ;
        
    % make sure P only contains unique points
    P = unique(P','rows')' ;
    
    % check if 3-D points are coplanar
    if size(P,1) == 3 && rank(P) < 3
        warning(['3-D points are coplanar! Projecting to first ',...
            'two dimensions!'])
        P = P(1:2,:) ;            
    end

    % STEP 3: take convex hull of points in workspace
    if size(P,2) > 1
        K = convhull(P') ;
    else
        K = 1 ;
    end
    
    % create output
    switch n_dim
        case 2
            P = P(:,K) ;
            n_P = size(P,2) ;
            F = [1:n_P,1] ;
        case 3
            F = K' ;
    end
    
    V = P' ;
end