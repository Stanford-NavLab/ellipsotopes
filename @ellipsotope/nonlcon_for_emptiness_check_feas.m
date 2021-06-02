function [c,ceq,gc,gceq] = nonlcon_for_emptiness_check_feas(E,x,p_norm,I)
% [c,ceq,gc,gceq] = nonlcon_for_emptiness_check_feas(E,x,p_norm,I)
%
% This function enforces the constraint that x lies in the ball product for
% the emptiness check corollary.
%
% Authors: Shreyas Kousik
% Created: 1 Jun 2021

    % get sizes of things
    n_x = length(x) ;
    n_I = length(I) ;
    
    % preallocate constraint and gradient
    c = zeros(n_I,1) ;
    gc = zeros(n_x,n_I) ;
    
    % compute constraint and gradient
    for idx = 1:n_I
        J = I{idx} ;
        c(idx) = sum(x(J).^p_norm) - 1 ;
        gc(J,idx) = p_norm*x(J).^(p_norm-1) ;
    end
    
    % equality constraint outputs
    ceq = [] ;
    gceq = [] ;
end