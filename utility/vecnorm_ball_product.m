function [V,V_mat] = vecnorm_ball_product(P_in,p_norm,I)
% V = vecnorm_ball_product(P_in,p_norm,I)
% [V,V_mat] = vecnorm_ball_product(P_in,p_norm,I)
%
% This function computes the generalized norm of points P_in \in \R^{n x
% m}, where n is the dimension and m is the number of points. The output V
% is given as
%
%   V(i) = max_{J \in I} ||P(J,i)||_p
%
% The optional second output, V_mat, is of the same size as P_in, with
%
%   V_mat(J,i) = ||P(J,i)||_p
%
% for each J \in I.
%
% The inputs are:
%   P_in     an n-by-m array of m n-dimensional points
%
%   p_norm   an even integer
%
%   I        an index set {J_1, J_2, ..., J_k} such that the union of all
%            J_i is the set {1,2,...,n} and, when i ~= j, the intersection of
%            J_i and J_j is empty.
%
% See also: project_points_to_ball_product
%
% Authors: Shreyas Kousik
% Created: 6 Apr 2021
% Updated: nup

% number of index subsets
n_I = length(I) ;

% initialize big output
V_mat = nan(size(P_in)) ;

% iterate through indices and do each projection
for idx_J = 1:n_I
    % get current index from the index set
    J = I{idx_J} ;
    
    % get the norm of all points for the current index set
    V_idx = vecnorm(P_in(J,:),p_norm,1) ;
    
    % assign the outputs
    V_mat(J,:) = repmat(V_idx,length(J),1) ;    
end

V = max(V_mat,[],1) ;