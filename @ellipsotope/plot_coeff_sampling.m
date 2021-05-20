function h_E = plot_coeff_sampling(E,varargin)
% plot_coeff_sampling(E)
% plot_coeff_sampling(E,'projdims',[dim1 dim2], other_input_args...)
% plot_coeff_sampling(E,'facecolor',color,'edgecolor',color,'facealpha',...)
%
% Plot an ellipsotope using coefficient sampling. This creates a patch
% object, and updates the E.plot_handle property.
%
% To plot an ellipsotope, we do the following:
%   1. generate a bunch of points in the coefficient space
%   2. push those points to the boundary of the feasible coefficient space
%   3. send the points through the affine map defined by the ellipsotope's
%      center and generators
%   4. plot the convex hull of the mapped points
%
% Authors: Adam Dai and Shreyas Kousik
% Created: in days of yore
% Updated: 27 Apr 2021 (updated emptiness warning message)

    %% prep/sanity check
    % get important properties
    [p,c,G,A,b,I] = E.get_properties() ;
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
    
    % if n_P is not set, set it to default value
    if exist('n_P','var') == 0
        n_P = min(7^size(G,2),10^5);
    end

    %% plotting setup
    % STEP 1: generate points in coefficient space
    % if the ellipsotope is basic...
    if E.is_basic()
        % switch how we generate points based on the coefficient space
        % dimension (using make_superllipse_2D, _3D, or _ND)
        if d_B == 2
            P = make_superellipse_2D(n_P,p);
        elseif d_B == 3
            P = make_superellipse_3D(n_P,p);
        else
            P = make_unit_superellipse_ND(n_P,p,d_B);
        end
    % else if...
    else
        % generate a bunch of random points 
        P = 2*rand(d_B,n_P) - 1 ;
        % project points to ball and linear subspace boundary (the existing
        % function will handle index sets and empty linear subspaces properly,
        % it turns out)
        [P,n_P] = project_points_to_ball_product_and_linear_subspace(P,p,A,b,I) ;
    end
    
    % check emptiness
    if isempty(P)
        warning('Ellipsotope to plot might be empty!')
    else
        % STEP 2: map points to ellipsotope workspace
        P = repmat(c,1,n_P) + G*P ;
        
        % make sure P only contains unique points
        P = unique(P','rows')' ;
        
        % STEP 3: take convex hull of points in workspace
        if size(P,2) > 1
            K = convhull(P') ;
        else
            K = 1 ;
        end
        
        n_K = length(K) ;
        
        %% plotting code here
        if n_K > 1
            plot_options = [{'facecolor','b','linewidth',1.0,'edgecolor','b',...
                'facealpha',0.1,'edgealpha',0.5}, varargin{:}] ;
        else
            plot_options = [{'b.','markersize',10}, varargin{:}] ;
        end
        
        switch length(proj_dims)
            case 1
                error('1-D plot is not implemented yet!')
            case 2
                P = P(:,K) ;
                P = points_to_CCW(P) ;
                n_P = size(P,2) ;
                if n_K > 1
                    h = patch('faces',[1:n_P,1],'vertices',P',plot_options{:}) ;
                else
                    h = plot(P(1),P(2),plot_options{:}) ;
                end
                
                E.plot_handle = h ;
            case 3
                if n_K > 1
                    h = trisurf(K,P(1,:)',P(2,:)',P(3,:)',plot_options{:}) ;
                else
                    h = plot(P(1),P(2),P(3),'b.','markersize',10) ;
                end
                
                E.plot_handle = h ;
        end
    end
    
    if nargout > 0
        h_E = h ;
    end
end