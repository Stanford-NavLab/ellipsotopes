function [out,value] = is_intersecting(E1,E2)
% [out,value] = is_intersecting(E1,E2)
%
% This function checks if two ellipsotopes are intersecting, by checking if
% ellipsotope formed from their intersection is empty. The output is 
% TRUE if the ellipsotopes are intersecting. The "value" output is the cost 
% of the emptiness convex program, which we set to 0 in the case
% of an unconstrained ellipsotope.
%
% Authors: Adam Dai
% Created: 9 May 2021
% Updated: 
    
    E_int = E1 & E2;
    [empty,value] = isempty(E_int);
    out = ~empty;
end