function plot(E,varargin)
% plot(E)
% plot(E,'projdims',[dim1 dim2], other_input_args...)
% plot(E,'facecolor',color,'edgecolor',color,'facealpha',...)
%
% Plot the ellipsotope if it is 2-D. This creates a patch
% object, and updates the E.plot_handle property.
%
% To plot an ellipsotope, we do the following:
%   1. generate a bunch of points in the coefficient space
%   2. push those points to the boundary of the feasible coefficient space
%   3. send the points through the affine map defined by the ellipsotope's
%      center and generators
%   4. plot the convex hull of the mapped points
%
% See also: test_plot_ellipsotope.m
%
% Authors: Adam Dai and Shreyas Kousik
% Created: in days of yore
% Updated: 20 Apr 2021 (minor cleanup)

    %% prep/sanity check
    % get important properties
    p = E.p_norm ;
    c = E.center ;
    G = E.generators ;
    A = E.constraint_A ; 
    b = E.constraint_b ; 
    I = E.index_set ;
    d = E.dimension ;
    d_B = size(G,2) ; % dimension of coefficient space

    % check for projection dimensions
    if nargin > 1 && strcmpi(varargin{1},'proj_dims')
        proj_dims = varargin{2} ;

        % remove from varargin
        if nargin > 3
            varargin = varargin{3:end} ;
        else
            varargin = {} ;
        end
    else
        % default to [1 2] projection
        proj_dims = 1:d ;
    end

    % check the dimension
    if d > 3 && ~exist('proj_dims','var')
        warning(['Plotting not supported for > 3-D ellipsotopes!',...
            'Plotting a 2-D projection instead'])
        % default to [1 2] projection
        proj_dims = 1:2 ;
    end

    % apply projection
    c = c(proj_dims) ;
    G = G(proj_dims,:) ;

    % check if E is basic, which allows us to reduce the
    % generator matrix nicely
    % if E.is_basic() && (p == 2) && size(G,1) <= size(G,2)
    %     % check if E is reduced
    %     if ~E.is_reduced()
    %         G = reduce_ellipsotope_generator_matrix(G) ;
    %     end
    % end

    % set the index set if it is empty
    if isempty(I)
        I = {1:d_B} ; % constrain ALL the coefficients!
    end

    %% plotting setup
    % STEP 1: generate points in coefficient space
    % if the ellipsotope is basic...
    if E.is_basic()
        % switch how we generate points based on the coefficient space
        % dimension (using make_superllipse_2D, _3D, or _ND)
        if d_B == 2
            n_P = 1000 ;
            P = make_superellipse_2D(n_P,p);
        elseif d_B == 3
            n_P = 1000 ;
            P = make_superellipse_3D(n_P,p);
        else
            n_P = 10000 ;
            P = make_unit_superellipse_ND(n_P,p,d_B);
        end
    % else if...
    else
        % generate a bunch of random points 
        n_P = 10000;
        P = 2*rand(d_B,n_P) - 1 ;
        % project points to ball and linear subspace boundary (the existing
        % function will handle index sets and empty linear subspaces properly,
        % it turns out)
        [P,n_P] = project_points_to_ball_product_and_linear_subspace(P,p,A,b,I) ;
    end
    
    % check emptiness
    if isempty(P)
        warning('Ellipsotope to plot is empty!')
    else
        % STEP 2: map points to ellipsotope workspace
        P = repmat(c,1,n_P) + G*P ;
        
        % STEP 3: take convex hull of points in workspace
        K = convhull(P') ;
        
        %% plotting code here
        patch_options = [{'facecolor','b','linewidth',1.5,'edgecolor','b',...
            'facealpha',0.1,'edgealpha',0.5}, varargin{:}] ;
        
        switch length(proj_dims)
            case 1
                error('1-D plot is not implemented yet!')
                % 2D
            case 2
                P = P(:,K) ;
                P = points_to_CCW(P) ;
                n_P = size(P,2) ;
                h = patch('faces',[1:n_P,1],'vertices',P',patch_options{:}) ;
                % 3D
            case 3
                h = trisurf(K,P(1,:)',P(2,:)',P(3,:)',patch_options{:}) ;
                
                E.plot_handle = h ;
        end
    end
end