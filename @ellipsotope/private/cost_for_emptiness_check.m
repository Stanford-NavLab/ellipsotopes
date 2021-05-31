function [c,gc] = cost_for_emptiness_check(E,x,p_norm,I,flag_bail_out_early)
% [c,gc] = E.cost_for_emptiness_check(x,p_norm,I,flag_bail_out_early)
%
% This is the (convex) cost function with subgradient for checking whether
% or not an ellipsotope is empty.
%
% Authors: Shreyas Kousik
% Created: 27 Apr 2021
% Updated: 30 May 2021 (minor speedups)

    % get number of index subsets
    n_I = length(I) ;
    
    % raise x to p
    x_pow = x.^p_norm ;
    
    % preallocate values to max over
    x_vals = nan(n_I,1) ;

    for idx = 1:n_I
        x_vals(idx) = sum(x_pow(I{idx})) ;
    end

    % compute cost
    [c,c_idx] = max(x_vals) ;
    
    if flag_bail_out_early && c <= 1
        error('Bailing out early!')
    end

    % compute cost (sub)gradient
    gc = zeros(1,length(x)) ;
    gc(I{c_idx}) = p_norm*(x(I{c_idx}).^(p_norm-1)) ;
end