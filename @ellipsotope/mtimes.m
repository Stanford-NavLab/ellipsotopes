% linear map (overloads *)
function out = mtimes(A, E)
% basic
if is_basic(E)
    out = ellipsotope(E.p_norm, A*E.center, A*E.generators);
    % constrained
elseif is_constrained(E)
    % generalized (and constrained)
    if is_general(E)
        out = ellipsotope(E.p_norm, A*E.center, A*E.generators, E.constraint_A, E.constraint_B, E.index_set);
        % not generalized (and constrained)
    else
        out = ellipsotope(E.p_norm, A*E.center, A*E.generators, E.constraint_A, E.constraint_B);
    end
    % not constrained
else
    % generalized (and not constrained)
    out = ellipsotope(E.p_norm, A*E.center, A*E.generators, [], [], E.index_set);
end
end