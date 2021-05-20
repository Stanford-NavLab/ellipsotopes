function h_E = plot(E,varargin)
% plot(E)
% plot(E,'projdims',[dim1 dim2], other_input_args...)
% plot(E,'facecolor',color,'edgecolor',color,'facealpha',...)
%
% Plot 2D and 3D ellipsotopes. This creates a patch object, and updates the 
% E.plot_handle property.
%
% Currently, there are two approaches for plotting:
%   1. plot_coeff_smapling.m
%       - Sampling in the coefficient space, and push these samples to the 
%         feasible boundary and affine through center and generators
%   2. plot_ray_tracing_2D.m
%       - Rotate a ray through angles to trace out the boundary of the
%         ellipsotope (only implemented for 2D)
% 
% The plot method selects between these two approaches based on the number
% of generators.
%
% See also: test_plot_ellipsotope.m
%
% Authors: Adam Dai and Shreyas Kousik
% Created: in days of yore
% Updated: 13 May 2021 (combined 2 plotting methods)

    %% prep/sanity check
    
    % default to coefficient method, unless in 2D and p==2
    d = E.dimension;
    p = E.p_norm;
    plot_method = 'coeff';
    if d == 2 && p == 2
        plot_method = 'ray';
    end
    
    %% TODO: move sanity checks here? 
    %%
    
    switch plot_method
        case 'coeff'
            h = plot_coeff_sampling(E,varargin{:});
        case 'ray'
            h = plot_ray_tracing_2D(E,varargin{:});
    end
    
    E.plot_handle = h;
    
    if nargout > 0
        h_E = h ;
    end

end