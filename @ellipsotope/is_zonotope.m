function out = is_zonotope(E)
% out = is_zonotope(E)
%
% Check if the ellipsotope E is a (constrained) zonotope.
%
% Authors: Shreyas Kousik
% Created: 20 May 2021
% Updated: 20 May 2021 (include inf-norm case)
    I = E.index_set ;
    G = E.generators ;
    n_gen = size(G,2) ;
    out = (length(I) == n_gen) || (E.p_norm == inf) ;
end