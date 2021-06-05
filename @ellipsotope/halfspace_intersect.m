function out = halfspace_intersect(E,h,f)
% halfspace_intersect(E,h,f)
% E.halfspace_intersect(h,f)
%
% Compute the intersection of an ellipsotopes with the halfspace defined by
%       h*x <= f
%
% Authors: Adam Dai 
% Created: 30 May 2021 
% Updated: 

% extract properties
[p,c,G,A,b,I] = E.get_properties;

% TODO: sanity check

% get sizes of things
n_gen = size(G,2) ;
n_con = size(A,1) ;

% perform intersection
d = 0.5*(sum(abs(h*G),2) + f - h*c) ;
G_S = [G, zeros(2,1)] ;
A_S = [A, zeros(n_con,1) ; h*G, d] ;
b_S = [b ; f - h*c - d] ;
I_S = [I, {n_gen + 1}] ;

out = ellipsotope(p,c,G_S,A_S,b_S,I_S);

end