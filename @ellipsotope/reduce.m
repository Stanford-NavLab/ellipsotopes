function out = reduce(E)
% reduce(E)
%
% Reduce the number of generators in E to dim(E.center).
    
    if E.is_basic()
        G_red = reduce_ellipsotope_generator_matrix(E.generators);
        out = ellipsotope(E.p_norm,E.center,G_red);
    else
        error('Order reduction is not yet implemented for non-basic ellipsotopes')
    end
end