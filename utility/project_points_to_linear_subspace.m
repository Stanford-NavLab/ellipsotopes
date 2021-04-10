function P_out = project_points_to_linear_subspace(P_in,A,b)
% P_out = project_points_to_linear_subspace(P_in,A,b)
%
% Given a linear subspace defined by Ax = b, and a bunch of points P_in,
% project the points onto the linear subspace.
%
% Authors: Shreyas Kousik
% Created: 9 Apr 2021
% Updated: nope

% get number of points
n_P = size(P_in,2) ;

% get null space and a feasible point
K = null(A) ;
t = A\b ;

% project points onto null space
P_dot_prods = K'*(P_in - repmat(t,1,n_P)) ;

% get the project points
P_out = K*P_dot_prods ;
end