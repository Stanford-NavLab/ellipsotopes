classdef ellipsotope < handle
    % class: ellipsotope
    %
    % This class implements the ellipsotope object
    %
    % Authors: Shreyas Kousik and Adam Dai
    % Created: 10 Feb 2021
    % Updated: 15 Mar 2021
    
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
        
        
        %% operations
        % linear map (overloads *)
        function out = mtimes(A, E)
            % basic
            if is_basic(E)
                out = ellipsotope(E.p_norm, A*E.center, A*E.generators);
            % constrained
            elseif is_constrained(E)
                % generalized (and constrained)
                if is_general(E)
                    out = ellipsotope(E.p_norm, A*E.center, A*E.generators, E.constraint_A, E.constraint_B, E.index_set);
                % not generalized (and constrained)
                else
                    out = ellipsotope(E.p_norm, A*E.center, A*E.generators, E.constraint_A, E.constraint_B);
                end
            % not constrained 
            else
                % generalized (and not constrained)
                out = ellipsotope(E.p_norm, A*E.center, A*E.generators, [], [], E.index_set);
            end
        end
        
        % intersection (overloads &)
        % assumes both inputs have same p_norm and are basic
        function out = and(E1, E2)
            c = E1.center;
            G = [E1.generators zeros(size(E1.generators))];
            A = [E1.generators -E2.generators];
            b = E2.center - E1.center;
            I = {1:E1.order,E1.order+1:E1.order+E2.order};
%             if is_constrained(E1)
%                 
%             else
%                 A = []; b= [];
%             end
%             if is_general(E)
%             else
%                 I = [];
%             end
            out = ellipsotope(E1.p_norm, c, G, A, b, I);
        end
        
        % minkowski addition (overloads +)
        % 
        function out = plus(s1, s2)
            % determine ellipsotope
            if isa(s1,'ellipsotope')
                out = s1;
                summand = s2;
            elseif isa(s2,'ellipsotope')
                out = s2;
                summand = s1;
            end
            
            % ellipsotope + ellipsotope
            if isa(summand,'ellipsotope')
                % both basic
                if is_basic(out) && is_basic(summand)
                    out.center = out.center + summand.center;
                    out.generators = [out.generators summand.generators];
                    out.index_set = {1:out.order,out.order+1:out.order+summand.order};
                % general case (constrained and generalized)
                else
                    out.center = out.center + summand.center;
                    out.generators = [out.generators summand.generators];
                    out.index_set = combine_indices(out.index_set, summand.index_set);
                    out.constraint_A = blkdiag(out.constraint_A, summand.constraint_A);
                    out.constraint_b = [out.constraint_b; summand.constraint_b];
                end
            % ellipsotope + vector
            elseif isnumeric(summand)
                out.center = out.center + summand;
            end
        end
        
        % containment
        %         function out = in(E, p)
        %             % out = in(E, p)
        %             %
        %             % Test if a point p \in \R^n is inside the ellipsotope E
        %         end
        
        %% plotting
        function plot(E,varargin)
            % plot(E)
            % plot(E,'projdims',[dim1 dim2], other_input_args...)
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
            
            % get important properties
            p = E.p_norm ;
            c = E.center ;
            G = E.generators ;
            d = E.dimension ;
            
            % set proj dims if dimension is greater than 3
            if ~exist('proj_dims','var')
                if d <= 3
                    proj_dims = 1:d ;
                else
                    proj_dims = 1:3 ;
                    warning(['Only plotting first ',num2str(d),' dimensions!'])
                end
            end
            G = G(proj_dims,:) ;
            
            % check if E is basic, which allows us to reduce the
            % generator matrix nicely
            if E.is_basic() && (p == 2)          
                % check if E is reduced
                if ~E.is_reduced()
                    G = reduce_ellipsotope_generator_matrix(G) ;
                end
            end
            
            % generate a bunch of points in the unit hypercube space and
            % plot the resulting object; we use "B" for "beta" which
            % we've used in our paper's notation for the ellipsotope
            % coefficients
            d_B = size(G,2) ; % dimension of coefficient unit hypercube
            switch d_B
                case 2
                    n_plot = 100 ;
                    B = make_superellipse_2D(p,1,zeros(2,1),n_plot) ;
                case 3
                    n_plot = 1000 ;
                    B = make_superellipse_3D(p,1,zeros(3,1),n_plot) ;
                otherwise
                    % get points on superellipse
                    n_plot = 10000 ; % we can be smarter than this I guess
                    B = make_unit_superellipse_ND(p,d_B,n_plot) ;
            end
            
            % check if we need to project the points onto a constraint
            if E.is_constrained()
                % get the constraints
                A = E.constraint_A ;
                b = E.constraint_b ;
                
                % get a basis for the constraint space
                N = null(A) ;
                b_offset = pinv(A)*b ;
                
                % project points onto the plane
                B = pinv(N*N')*(B - b_offset) + b_offset ;
                
                % add some more points in the nullspace (shrugs)
                n_extra = 10000 ;
                B_extra = 4*rand(size(N,2),n_extra) - 2 ;
                B = [B, N*B_extra + b_offset] ;
                
                % keep the points that obey the norm
                for idx = 1:length(E.index_set)
                    N_log = vecnorm(B(E.index_set{idx},:),p) <= 1 ;
                    B = B(:,N_log) ;
                end
            end
            
            % map points to the ellipsotope's proj dims (recall that the
            % generators have already been projected)
            P = G*B + repmat(c,1,size(B,2)) ;
            
            % get the convex hull of the points
            ch = convhull(P') ;
            
            % set default plot inputs
            patch_options = [{'facecolor','b','edgecolor','b',...
                    'facealpha',0.1,'edgealpha',1}, varargin] ;
            
            % plot!
            switch length(proj_dims)
                case 1
                    error('1-D plot is not implemented yet!')
                case 2
                    P = P(:,ch) ;
                    P = points_to_CCW(P) ;
                    F = [1:size(P,2),1] ;
                    h = patch('faces',F,'vertices',P',patch_options{:}) ;
                case 3
                    h = trisurf(ch,P(1,:)',P(2,:)',P(3,:)',patch_options{:}) ;
            end
            
            % finish up
            E.plot_handle = h ;
        end
    end
end