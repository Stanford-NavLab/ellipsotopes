function E_rdc = reduce_constrained_2_etope(E)
% E_rdc = reduce_constrained_2_etope(E)
%
% Turn a constrained 2-ellipsotope into a basic 2-ellipsotope by leveraging
% the fact that they are both ellipsoids
%
% Authors: Shreyas Kousik
% Created: 13 July 2021
% Updated: 14 July 2021 (took this out of the @ellipsotope methods)

    % get relevant properties
    [p,c,G,A,b,~,~,~,~,n_I] = get_properties(E) ;

    if p ~= 2
        error('This function only works for 2-ellipsotopes!')
    end

    if ~isempty(A)
        if n_I > 1
            error(['Cannot use this method for an indexed 2-ellipsotope! ',...
                'Try E.identify_component_ellipsotopes first!'])
        else
            % get affine map to hyperplane
            t = pinv(A)*b ; % center of intersection area
            rd = sqrt(1 - vecnorm(t)^2) ; % radius of intersected ball
            K = null(A) ;
            T = rd*K ;

            % compute new etope matrix
            c_rdc = c + G*t ;
            G_rdc = G*T ;
        end
    else
        c_rdc = c;
        G_rdc = G ;
    end

    % get minimal G_rdc
    G_rdc = reduce_2_etope_generator_matrix(G_rdc) ;

    % create new basic etope
    E_rdc = ellipsotope(2,c_rdc,G_rdc) ;
end