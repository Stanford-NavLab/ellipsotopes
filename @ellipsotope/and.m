function out = and(E1, E2)
% and(E1, E2)
% E1 & E2
%
% Compute the intersection of 2 ellipsotopes (as a new ellipsotope). 
% Overloads the & operator for 2 ellipsotopes.
% Assumes both ellipsotopes have the same p_norm and are basic.
%
% Authors: Adam Dai 
% Created: 1 Mar 2021 
% Updated: 

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