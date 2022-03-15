function [F,V] = make_patch_data_constrained_zonotope(c,G,A,b)
% [F,V] = make_patch_data_constrained_zonotope(c,G,A,b)
%
% Given a constrained zonotope Z = (c,G,A,b), return a set of faces and
% vertices so that it can be plotted in either 2-D or 3-D.
%
% Added a sanity check for large n_gen!
%
% Authors: Shreyas Kousik
% Created: 24 Mar 2021
% Updated: 15 Mar 2022

    % get dimemsion
    n_gen = size(G,2) ;
    
    if n_gen > 17
        error('This is a lot of generators! You probably don''t want to try to plot this zonotope!')
    end

    % construct inequality constraints of hypercube
    A_ineq = [eye(n_gen) ; -eye(n_gen)] ;
    b_ineq = ones(2*n_gen,1) ;

    % try to get vertices
    try
        V_idx = lcon2vert(A_ineq,b_ineq,A,b) ;
        V_idx = c + G*V_idx' ;
        try
            F_idx = convhull(V_idx') ;
        catch
            F_idx = [] ;
        end
        V = V_idx' ;
        F = F_idx' ;
    catch
        % set default output
        F = [] ;
        V = [] ;
    end
end