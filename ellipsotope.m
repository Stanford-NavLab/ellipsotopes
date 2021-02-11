classdef ellipsotope < handle
% class: ellipsotope
%
% This class implements the ellipsotope object
%
% Authors: Shreyas Kousik and Adam Dai
% Created: 10 Feb 2021
    properties
        % basic ellipsotope properties
        p_norm_ = 2 ; % \in 2*N, default is p = 2, zonotope is p = Inf
        center = [] ; % c \in \R^n
        generators = [] ; % [g_1, g_2, ..., g_m] \in \R^(n x m)
        
        % constrained ellipsotope properties
        constraint_A = [] ; % A*beta = b, where beta is the generator coeff
        constraint_b = [] ; 
        
        % generalized ellipsotope
        indices = {} ; % e.g., if m = 3, then I = {[1,2],[2,3]}
        
        
        order = [] ; % number of generators
    end
    
    
    %% methods
    methods
        %% constructor
        function E = ellipsotope(p,c,G,A,b,I)
            
        end
    end
end