function [p,c,G,A,b,I] = get_properties(E)
% [c,G,A,b,I] = get_properties(E)
%
% Get all the defining properties of the ellipsotope E
%
% Authors: Shreyas Kousik
% Created: 21 Apr 2021
% Updated: nah bruh
p = E.p_norm ;
c = E.center ;
G = E.generators ;
A = E.constraint_A ;
b = E.constraint_b ;
I = E.index_set ;
end