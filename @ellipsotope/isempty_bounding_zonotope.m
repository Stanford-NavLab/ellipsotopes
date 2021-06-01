function [out,value] = isempty_bounding_zonotope(E)
% [out,value] = isempty_bounding_zonotope(E)
%
% Run an LP with the constrained zonotope that outer-bounds an ellipsotope
% to check if it is empty conservatively
%
% The output is TRUE if the ellipsotope is EMPTY, in which case value <= 1.
%
% Authors: Shreyas Kousik
% Created: 1 June 2021
% Updated: nah
    [f_cost,A_ineq,b_ineq,A_eq,b_eq] = E.make_con_zono_empty_check_LP(E.constraint_A,E.constraint_b) ;
    z_opt = linprog(f_cost,A_ineq,b_ineq,A_eq,b_eq)  ;
    value = z_opt(end) ;
    out = value > 1 ;
end