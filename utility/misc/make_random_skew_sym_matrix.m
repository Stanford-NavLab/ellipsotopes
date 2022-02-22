function S = make_random_skew_sym_matrix(n)
% S = make_random_skew_sym_matrix(n)
%
% Makes an n-dimensional random skew-symmetric matrix.
%
% Authors: Shreyas Kousik
% Created: 21 Feb 2022
% Updated: nah
    S = 2*rand(n) - 1 ;
    S = triu(S) ;
    S = S - S' ;
end