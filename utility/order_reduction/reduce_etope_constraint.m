function [E_rdc,idx_con_rdc] = reduce_etope_constraint(E,idx_con_rdc)
% E = reduce_etope_constraint(E,idx_rdc)
%
% Just reduce a single constraint (no generator) for an ellipsotope (this
% works for any p-norm).
%
% Authors: Shreyas Kousik
% Created: 16 Mar 2022
% Updated: --

    [p,c,G,A,b,I,~,~,n_con] = get_properties(E) ;

    if nargin < 2
        % get the generator associated with the best reduction
        idx_rdc = estimate_best_idx_for_gen_and_con_removal(E) ;
        
        % find a constraint where the j-th entry is nonzero
        a_ij_all = A(:,idx_rdc) ;
        idx_con_rdc = find(a_ij_all ~= 0,1) ;
    end

    % create Lambda
    lm = zeros(n_con,1) ;
    lm(idx_con_rdc) = 1 ;
    Lm = diag(lm) ;

    % create new etope
    A_rdc = A - Lm*A ;
    b_rdc = b - Lm*b ;

    % delete constraint
    A_rdc = A_rdc(~lm,:) ;
    b_rdc = b_rdc(~lm) ;

    E_rdc = ellipsotope(p,c,G,A_rdc,b_rdc,I) ;
end