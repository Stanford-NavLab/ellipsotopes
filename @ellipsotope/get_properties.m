function [p,c,G,A,b,I,n_dim,n_gen,n_con,n_I] = get_properties(E)
% [p,c,G,A,b,I] = E.get_properties()
% [p,c,G,A,b,I,n_dim,n_gen,n_con] = E.get_properties()
%
% Get all the defining properties of the ellipsotope E
%
% Authors: Shreyas Kousik
% Created: 21 Apr 2021
% Updated: 21 Sep 2021 (edited help text)

p = E.p_norm ;
c = E.center ;
G = E.generators ;
A = E.constraint_A ;
b = E.constraint_b ;
I = E.index_set ;

[n_dim,n_gen] = size(G) ;
n_con = size(A,1) ;

n_I = length(I) ;
end