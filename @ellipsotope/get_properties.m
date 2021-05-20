function [p,c,G,A,b,I,n_dim,n_gen,n_con] = get_properties(E)
% [c,G,A,b,I] = E.get_properties()
% [p,c,G,A,b,I,n_dim,n_gen,n_con] = E.get_properties()
%
% Get all the defining properties of the ellipsotope E
%
% Authors: Shreyas Kousik
% Created: 21 Apr 2021
% Updated: 20 May 2021 (edited help text, added more output args)
p = E.p_norm ;
c = E.center ;
G = E.generators ;
A = E.constraint_A ;
b = E.constraint_b ;
I = E.index_set ;

[n_dim,n_gen] = size(G) ;
n_con = size(A,1) ;

end