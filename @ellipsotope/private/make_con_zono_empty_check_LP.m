function [f_cost,A_ineq,b_ineq,A_eq,b_eq] = make_con_zono_empty_check_LP(E,A,b)
% [f_cost,A_ineq,b_ineq,A_eq,b_eq] = E.make_con_zono_empty_check_LP(A,b)
%
% Given the constraint matrices A and b for a constrained zonotope, return
% the data matrices required to construct a linear program to solve the
% program in [1, Eq. (19)]:
%
%   min{ max{|x|} : Ax = b}
%
% References:
% [1] Constrained zonotopes, https://web.mit.edu/braatzgroup/Scott_Automatica_2016.pdf
%
% Authors: Shreyas Kousik
% Created: 30 May 2021

    % dimension of problem
    d = size(A,2) ;

    % cost
    f_cost = [zeros(d,1); 1] ;

    % inequality cons
    A_ineq = [-eye(d), -ones(d,1) ; eye(d), -ones(d,1)] ;
    b_ineq = zeros(2*d,1) ;

    % equality cons
    A_eq = [A, zeros(size(A,1),1)] ;
    b_eq = b ;
end