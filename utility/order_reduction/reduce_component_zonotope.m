function E_rdc = reduce_component_zonotope(E,n_rdc)
% E_rdc = reduce_component_zonotope(E,n_rdc)
%
% Identify and reduce component zonotopes by lifting E and reducing. Works
% for any p!
%
% Authors: Shreyas Kousik
% Created: 16 Mar 2022
% Updated: nah

    % lift!
    n_dim = E.dimension ;
    E_l_orig = lift(E) ;
    [E_l,idx_Z_I,idx_Z] = identify_component_zonotope(E_l_orig) ;
    [p,c_l,G_l,~,~,I_l,~,n_gen,~] = get_properties(E_l) ;

    % figure out how much to reduce the zonotope generators
    if nargin < 2
        n_rdc = n_gen - n_des ;
    end

    % get the lifted generator matrix of the component zonotope
    G_rdc = G_l(:,idx_Z:n_gen) ;
    [n_r,n_c] = size(G_rdc) ;

    % if G_rdc is "wide" then it can be reduced, otherwise we need to
    % pop some generators...

    if n_c > n_r
        G_rdc = reduce_zonotope_Chischi(G_rdc,n_rdc) ;

        % create a new index set for the reduced generator matrix
        I_rdc = I_l(1:(idx_Z_I-1)) ;
        idx_max = get_max_index(I_rdc) ;
        I_rdc = [I_rdc, num2cell((1:size(G_rdc,2)) + idx_max)] ;

        % reconstruct tope
        E_l = ellipsotope(p,c_l,[G_l(:,1:(idx_Z-1)),G_rdc],[],[],I_rdc) ;
        E_rdc = drop(E_l,n_dim) ;
    end
end