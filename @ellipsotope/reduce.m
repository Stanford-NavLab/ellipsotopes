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