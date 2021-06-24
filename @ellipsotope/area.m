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
        % call plot to store vertices in plot_handle
        plot(E,'plot_method','ray');

        % retrieve vertices, and use to compute area
        V = E.plot_handle.Vertices;
        out = area(alphaShape(V));
    end
end