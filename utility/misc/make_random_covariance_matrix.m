function S = make_random_covariance_matrix(d)
% S = make_random_covariance_matrix(d)
%
% It does what it says on the box. The input d is the dimension, so the
% output is a d x d matrix.
%
% Authors: Shreyas Kousik
    Q = rand(d);
    D = diag(rand(d,1)) ;
    S = Q*D*Q' ;
end