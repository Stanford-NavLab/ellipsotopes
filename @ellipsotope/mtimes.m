function out = mtimes(A, E)
% mtimes(A, E)
% A * E
%
% Compute linear mapping of ellipsotope under a matrix.
% Overloads the * operator for 1 ellipsotope and 1 matrix.
%
% Authors: Adam Dai 
% Created: 1 Mar 2021 
% Updated: 25 Mar 2021

% basic
if is_basic(E)
    out = ellipsotope(E.p_norm, A*E.center, A*E.generators);
% constrained
elseif is_constrained(E)
    % indexed (and constrained)
    if is_general(E)
        out = ellipsotope(E.p_norm, A*E.center, A*E.generators, E.constraint_A, E.constraint_B, E.index_set);
        % not indexed (and constrained)
    else
        out = ellipsotope(E.p_norm, A*E.center, A*E.generators, E.constraint_A, E.constraint_B);
    end
    % not constrained
else
    % indexed (and not constrained)
    out = ellipsotope(E.p_norm, A*E.center, A*E.generators, [], [], E.index_set);
end
end