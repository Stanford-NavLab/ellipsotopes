function out = is_constrained(E)
out = ~((isempty(E.constraint_A)) || ...
    (isempty(E.constraint_b)));
end