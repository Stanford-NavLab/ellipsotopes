function [c,gc] = cost_for_emptiness_check_feas(E,x,A,b)
% [c,gc] = E.cost_for_emptiness_check_feas(x,A,b)
%
% This cost function attempts to find a feasible point in an ellipsotope
% using the corollary in the paper.
%
% Authors: Shreyas Kousik
% Created: 1 Jun 2021
    c = sum((A*x - b) .^2) ;
    gc = (2*(A'*A)*x)' - 2*b'*A ;
end
