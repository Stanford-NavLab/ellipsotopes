function h_E = plot(E,varargin)
% plot(E)
% plot(E,'projdims',[dim1 dim2], other_input_args...)
% plot(E,'plot_method',method, other_input_args...)
% plot(E,'facecolor',color,'edgecolor',color,'facealpha',...)
%
% Plot 2D and 3D ellipsotopes. This creates a patch object, and updates the 
% E.plot_handle property.
%
% In addition to the usual patch keywork input arguments, one can pass in
% 'proj_dims' to pick which projected dimensions of the ellipsotope to
% plot.
%
% One can also pass in 'plot_method', followed by a string:
%
%   'sample' - this samples in the ellipsotope coefficient space, pushes
%              the samples to the feasible boundary in that space, then
%              maps them to the ellipsotope workspace
%
%   'ray'    - this rotates a ray in 2-D through 2*pi radians to trace out
%              the ellipsotope boundary
%
%   'zono'   - this is used by default if the ellipsotope is a zonotope, or
%              plots the constrained zonotope enclosing an ellipsotope if
%              specified by the user
%
% If no plot method is specified, we choose which one to use based on the
% ellipsotope dimension and number of generators.
%
% Finally, one can pass in 'num_points', which specifies how many points to
% plot using either the sample or ray methods. We recommend not specifying
% this for the sample method, which sets the number of points dynamically
% based on the number of generators. For the ray tracing method, the
% default is 200.
%
% See also: test_plot_ellipsotope.m
%
% Authors: Adam Dai and Shreyas Kousik
% Created: in days of yore
% Updated: 20 May 2021 (restructured all plotting)

    %% setup
    % get properties
    [p,c,G,A,b,I,n_dim,n_gen,~] = E.get_properties() ;
    
    % check for projection dimensions
    [proj_dims,args] = check_varargin_for_keyword('proj_dims',varargin{:}) ;
    if isempty(proj_dims)
        % default to [1 2] projection
        proj_dims = 1:2 ;
    end
    % project
    c = c(proj_dims) ;
    G = G(proj_dims,:) ;
    
    % check for user-specified method
    [plot_method,args] = check_varargin_for_keyword('plot_method',args{:}) ;
    if isempty(plot_method)
        % default to coefficient method, unless in 2D
        plot_method = 'sample';
        if (n_dim == 2) && (p == 2)
            plot_method = 'ray';
        end
    end
    
    % check for number of points
    [n_P,args] = check_varargin_for_keyword('num_points',args{:}) ;
    
    % create patch data args in
    patch_data_args_in = {p,c,G,A,b,I} ;
    if ~isempty(n_P) 
        patch_data_args_in = [patch_data_args_in,{n_P}] ;
    end
    
    % override plot method if E is a zonotope
%     if E.is_zonotope()
%         plot_method = 'zono' ;
%         patch_data_args_in = {c,G,A,b} ;
%     end    
    
    %% create plot data
    switch plot_method
        case 'sample'           
            [F,V] = make_patch_data_coeff_sampling(patch_data_args_in{:}) ;
        case 'ray'
            [F,V] = make_patch_data_ray_tracing_2D(patch_data_args_in{:}) ;
        case 'zono'
            [F,V] = make_patch_data_constrained_zonotope(patch_data_args_in{:}) ;
        otherwise
            error('Invalid plot method! Please choose sample, ray, or zono.')
    end
    
    %% plotting
    plot_options = [{'facecolor','b','linewidth',1.5,'edgecolor','b',...
                'facealpha',0.1,'edgealpha',0.5}, args{:}] ;
    h = patch('faces',F,'vertices',V,plot_options{:}) ;
    
    E.plot_handle = h;
    
    if nargout > 0
        h_E = h ;
    end

end