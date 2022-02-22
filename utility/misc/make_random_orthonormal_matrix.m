function [R,S] = make_random_orthonormal_matrix(n)
% R = make_random_orthonormal_matrix(n)
%
% Authors: Shreyas Kousik
% Created: 21 Feb 2022
% Updated: nah

    % create a random skew-symmetrix matrix
    S = pi.*make_random_skew_sym_matrix(n) ;

    % exponentiate!
    R = expm(S) ;
end