classdef ellipsotope < handle
    % class: ellipsotope
    %
    % This class implements the ellipsotope object
    %
    % Authors: Shreyas Kousik and Adam Dai
    % Created: 10 Feb 2021
    properties
        % basic ellipsotope properties
        p_norm = 2 ; % \in 2*N, default is p = 2, zonotope is p = Inf
        center = [] ; % c \in \R^n
        generators = [] ; % [g_1, g_2, ..., g_m] \in \R^(n x m)
        ellipsoid_shape = [] ;
        
        % constrained ellipsotope properties
        constraint_A = [] ; % A*beta = b, where beta is the generator coeff
        constraint_b = [] ;
        
        % generalized ellipsotope
        index_set = {} ; % e.g., if m = 3, then I = {[1,2],[2,3]}
        
        % complexity
        dimension = [] ;
        order = [] ;
        
        % plot methods
        plot_handle = [] ;
    end
    
    
    %% methods
    methods
        %% constructor
        function E = ellipsotope(p,c,G,A,b,I)
            % E = ellipsotope(p,c,G)
            %
            % Create an ellipsotope. By default, this just assumes a basic
            % ellipsotope.
            %
            % To create a constrained ellipsotope, use
            %   E = ellipsotope(p,c,G,A,b)
            %
            % To create a general ellipsotope, use
            %   E = ellipsotope(p,c,G,[],[],I)
            % where I is a set of indices
            %
            % Authors: Shreyas Kousik
            % Created: 10 Feb 2021
            
            % basic ellipsotope properties
            E.p_norm = p ;
            E.center = c ;
            E.generators = G ;
            
            if nargin > 3
                E.constraint_A = A ;
                E.constraint_b = b ;
            end
            
            if nargin > 5
                E.index_set = I ;
            end
            
            % get the size
            E.dimension = size(G,1) ;
            E.order = size(G,2) ;
            
            % get the ellipsoid shape matrix
            if E.is_reduced()
                
            end
        end
        
        %% property tests
        function out = is_basic(E)
            % out = is_basic(E)
            %
            % Test if E is a basic ellipsotope, and return true/false.
            out = ~(E.is_constrained || E.is_general) ;
        end
        
        function out = is_constrained(E)
            out = ~((isempty(E.constraint_A)) || ...
                (isempty(E.constraint_b)));
        end
        
        function out = is_general(E)
            out = ~isempty(E.index_set) ;
        end
        
        function out = is_reduced(E)
            d = E.dimension ;
            o = E.order ;
            
            out = E.is_basic && (d == o) ;
        end
        
        %% order reduction
        function reduce(E)
            % reduce(E)
            %
            % Reduce the number of generators in E to dim(E.center).
            if E.is_basic()
                E.generators = reduce_ellipsotope_generator_matrix(E.generators) ;
                E.order = size(E.generators,2) ;
            else
                error('Order reduction is not yet implemented for non-basic ellipsotopes')
            end
        end
        
        %% plotting
        function plot(E,varargin)
            % plot(E)
            % plot(E,'projdims',[dim1 dim1], other_input_args...)
            % plot(E,'facecolor',color,'edgecolor',color,'facealpha',...)
            %
            % Plot the ellipsotope if it is 2-D. This creates a patch
            % object, and updates the E.plot_handle property.
            
            % check the dimension
            if E.dimension > 2
                warning(['Plotting not supported for > 2-D ellipsotopes!',...
                    'Plotting a 2-D projection instead'])
            end
            
            % check for projection dimensions
            if nargin > 1 && strcmpi(varargin{1},'projdims')
                proj_dims = varargin{2} ;
                
                % clean up varargin
                if length(varargin) > 2
                    varargin = varargin(3,:) ;
                end
            end
            
            % check if E is basic
            if E.is_basic()
                G = E.generators ;
                
                % check if E is reduced
                if E.is_reduced()
                    G = reduce_ellipsotope_generator_matrix(G) ;
                end
                
                % make a circle of points
                [F,V] = make_circle(1,100) ;
                
                % map the vertices using the generator matrix (amazing!)
                V = (G*V')' ;
                
                % plot the ellipsotope
                if nargin == 1
                    h = patch('faces',F,'vertices',V,...
                        'facecolor','b','edgecolor','b',...
                        'facealpha',0.1,'edgealpha',1) ;
                else
                    h = patch('faces',F,'vertices',V,varargin{:}) ;
                end
                
                E.plot_handle = h ;
            else
                error('Plotting for non-basic ellipsotopes is not working yet')
            end
        end
    end
end