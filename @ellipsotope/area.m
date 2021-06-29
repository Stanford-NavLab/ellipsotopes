function out = area(E)
% out = area(E)
%
% Computes approximate area of a 2D ellipsotope
    
    % check that e-tope is 2D
    assert(E.dimension == 2, 'Cannot compute area of non-2D Ellipsotope')
    
    % if zonotope
    if E.is_zonotope()
        cz = conZonotope(E.center,E.generators);
        V = vertices(cz);
        out = area(alphaShape(V(1,:)',V(2,:)'));
    else
        % using ray tracing to obtain "vertices" (boundary points)
        n_P = 200;
        [p,c,G,A,b,I,~,~,~] = get_properties(E);
        [~,V] = make_patch_data_ray_tracing_2D(p,c,G,A,b,I,n_P);

        % create polyshape from vertices and use to compute area
        out = area(polyshape(V));
    end
end