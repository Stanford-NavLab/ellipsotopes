function h_E = plot(E,varargin)
% plot(E)
% plot(E,'proj_dims',[dim1 dim2], other_input_args...)
% plot(E,'plot_method',method, other_input_args...)
% plot(E,'facecolor',color,'edgecolor',color,'facealpha',...)
%
% Plot 2D and 3D ellipsotopes. This creates a patch object, and updates the
% E.plot_handle property.
%
% In addition to the usual patch keywork input arguments, one can pass in
% the following special keywords:
%   'proj_dims' - pick which projected dimensions to plot
%
%   'plot_all_dims' - create as many 2-D subplots as needed to plot all
%                     projected dims of the ellipsotope
%
%   'color' - set the facecolor and edgecolor arguments in one shot
%
%   'plot_method' - see below
%
%   'num_points' - how many points to plot using sample or ray tracing
%                  methods; the sample method sets the number dynamically,
%                  whereas ray tracing defaults to 200
%
% If you pass in 'plot_method', also pass in a string from the following:
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
% See also: test_plot_ellipsotope.m
%
% Authors: Adam Dai and Shreyas Kousik
% Created: in days of yore
% Updated: 14 Mar 2022 (added 'plot_all_dims')

%% setup
% check emptiness
if E.isempty()
    warning('Ellipsotope is empty! Not plotting!')
else
    % get properties
    [p,c,G,A,b,I,n_dim,~,~] = E.get_properties() ;
    
    % check for projection dimensions
    [flag_plot_all_dims,plot_args] = check_varargin_for_keyword('plot_all_dims',varargin{:}) ;
    
    if ~isempty(flag_plot_all_dims) && flag_plot_all_dims
        % plot each pair of dimensions
        idx_dims = 1:2:n_dim ; % all dimensions
        n_plots = length(idx_dims) ;
        idx_plot = 1 ;
        
        h = [] ;
        for idx = idx_dims
            % make sub-tope to plot
            if idx == n_dim
                proj_dims = [idx-1, idx] ;
            else
                proj_dims = [idx, idx+1] ;
            end
            c_idx = c(proj_dims) ;
            G_idx = G(proj_dims,:) ;
            
            E_idx = ellipsotope(p,c_idx,G_idx,A,b,I) ;
            
            subplot(1,n_plots,idx_plot) ; grid on ; axis equal ;
            args_in = [plot_args, {'plot_all_dims',false}] ;
            h_idx = plot(E_idx,args_in{:}) ;
            xlabel(['dim ',num2str(proj_dims(1))])
            ylabel(['dim ',num2str(proj_dims(2))])
            h = [h,h_idx] ;
            idx_plot = idx_plot + 1 ;
        end
    else
        [proj_dims,plot_args] = check_varargin_for_keyword('proj_dims',plot_args{:}) ;
        if isempty(proj_dims)
            % default to [1 2] projection
            proj_dims = 1:2 ;
            
            if n_dim == 3
                proj_dims = 1:3 ;
            end
        end
        
        % project
        c = c(proj_dims) ;
        G = G(proj_dims,:) ;
        
        % check for user-specified method
        [plot_method,plot_args] = check_varargin_for_keyword('plot_method',plot_args{:}) ;
        if isempty(plot_method)
            % default to coefficient method, unless in 2D
            plot_method = 'sample';
            if (n_dim == 2) && (p == 2)
                plot_method = 'ray';
            end
        end
        
        % check for color input argument
        [color,plot_args] = check_varargin_for_keyword('color',plot_args{:}) ;
        if ~isempty(color)
            plot_args = [plot_args,...
                {'facecolor',color,...
                'edgecolor',color}] ;
        end
        
        % check for number of points
        [n_P,plot_args] = check_varargin_for_keyword('num_points',plot_args{:}) ;
        
        % create patch data args in
        patch_data_args_in = {p,c,G,A,b,I} ;
        if ~isempty(n_P)
            patch_data_args_in = [patch_data_args_in,{n_P}] ;
        end
        
        % override plot method if E is a zonotope
        if E.is_zonotope()
            plot_method = 'zono' ;
            patch_data_args_in = {c,G,A,b} ;
        end
        
        %% create plot data
        switch plot_method
            case 'sample'
                [F,V] = make_patch_data_coeff_sampling(patch_data_args_in{:}) ;
            case 'ray'
                [F,V] = make_patch_data_ray_tracing_2D(patch_data_args_in{:}) ;
            case 'zono'
                try
                    [F,V] = make_patch_data_constrained_zonotope(patch_data_args_in{:}) ;
                catch
                    warning(['Cannot plot a zonotope with lots of generators! ',...
                        'Using ray tracing method instead!'])
                    % create patch data args in
                    patch_data_args_in = {p,c,G,A,b,I} ;
                    if ~isempty(n_P)
                        patch_data_args_in = [patch_data_args_in,{n_P}] ;
                    end
                    [F,V] = make_patch_data_ray_tracing_2D(patch_data_args_in{:}) ;
                end
            otherwise
                error('Invalid plot method! Please choose sample, ray, or zono.')
        end
        
        %% plotting
        plot_options = [{'facecolor','b','linewidth',1.5,'edgecolor','b',...
            'facealpha',0.1,'edgealpha',0.5}, plot_args{:}] ;
        h = patch('faces',F,'vertices',V,plot_options{:}) ;
        
        E.plot_handle = h;
        
        if nargout > 0
            h_E = h ;
        end
    end
end
end