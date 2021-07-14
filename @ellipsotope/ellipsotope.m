classdef ellipsotope < handle
    % class: ellipsotope
    %
    % This class implements the ellipsotope object
    %
    % Authors: Shreyas Kousik and Adam Dai
    % Created: 10 Feb 2021
    % Updated: 14 July 2021 (added constraint cleanup)
    
    properties
        % basic ellipsotope properties
        p_norm = 2 ; % \in 2*N, default is p = 2, zonotope is p = Inf
        center = [] ; % c \in \R^n
        generators = [] ; % [g_1, g_2, ..., g_m] \in \R^(n x m)
        ellipsoid_shape = [] ;
        
        % constrained ellipsotope properties
        constraint_A = [] ; % A*beta = b, where beta is the generator coeff
        constraint_b = [] ;
        
        % indexed ellipsotope
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
            % To create an indexed ellipsotope, use
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
                % remove all zero rows
                log_A_zero_rows = all(A == 0,2) ;
                A = A(~log_A_zero_rows,:) ;
                b = b(~log_A_zero_rows) ;
                
                % set constraints
                E.constraint_A = A ;
                E.constraint_b = b ;
            end
            
            if nargin > 5
                if ~check_index_set_validity(I,G)
                    error('The index set is invalid!')
                end
                E.index_set = I ;
            end
            
            % get the size
            E.dimension = size(G,1) ;
            E.order = size(G,2) ;
            
            % if index set is empty, initialize with default index set
            if isempty(E.index_set)
                E.index_set = {1:E.order};
            end
            
            % get the ellipsoid shape matrix (TODO: unfinished?)
            if E.is_reduced()
                
            end
        end
    end
end