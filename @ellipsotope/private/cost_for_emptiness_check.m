function [c,gc] = cost_for_emptiness_check(E,x,p_norm,I)
% [c,gc] = E.cost_for_emptiness_check(x,p_norm,I)
%
% This is the (convex) cost function with subgradient for checking whether
% or not an ellipsotope is empty.
%
% Authors: Shreyas Kousik
% Created: 27 Apr 2021
% Updated: nope

    % get number of index subsets
    n_I = length(I) ;
    
    % preallocate values to max over
    x_vals = nan(n_I,1) ;

    for idx = 1:n_I
        x_vals(idx) = sum(x(I{idx}).^p_norm) ;
    end

    % compute cost
    [c,c_idx] = max(x_vals) ;

    % compute cost (sub)gradient
    gc = zeros(1,length(x)) ;
    gc(I{c_idx}) = p_norm*(x(I{c_idx}).^(p_norm-1)) ;
end