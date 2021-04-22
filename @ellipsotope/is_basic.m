function out = is_basic(E)
% out = is_basic(E)
%
% Test if E is a basic ellipsotope, and return true/false.
out = ~(E.is_constrained || E.is_indexed) ;
end