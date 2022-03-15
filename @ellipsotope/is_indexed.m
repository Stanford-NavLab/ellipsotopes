function out = is_indexed(E)
% true_or_false = is_indexed(E)
%
% Authors: Shreyas
% Created: shrug
% Updated: 14 Mar 2022 (fixed to check for constraints)

    out = (length(E.index_set) > 1) && (~E.is_constrained()) ;
end