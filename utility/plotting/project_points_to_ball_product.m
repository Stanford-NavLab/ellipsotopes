function P_out = project_points_to_ball_product(P_in,p_norm,I,idxs_J_to_enf)
% P_out = project_points_to_ball_product(P_in,p_norm,I)
% P_out = project_points_to_ball_product(P_in,p_norm,I,idxs_J_to_enf)
%
% This function projects the input points P_in to the boundary of the ball
% product defined by p_norm and the index set I. The optional input
% idx_J_to_enf picks which multiindex in the index set to enforce for each
% point; it is created at random if not provided.
%
% See also: test_projection_to_ball_product,
%           project_points_to_ball_product_and_linear_subspace
%
% Authors: Shreyas Kousik
% Created: 5 Apr 2021
% Updated: 13 Apr 2021

% sizes of things
n_I = length(I) ;
n_P = size(P_in,2) ;

% set default args in
if nargin < 4
    % make random list of which norm to enforce per point
    idxs_J_to_enf = rand_int(1,n_I,[],[],1,n_P) ;
end

% initialize output
P_out = P_in ;

% iterate through indices and do each projection
for idx_J = 1:n_I
    % get current index from the index set
    J = I{idx_J} ;
    
    % get the norm of all points for the current index set
    V = vecnorm(P_out(J,:),p_norm,1) ;
    
    % get all points on which to enforce the norm for this index, which
    % includes all points for which the norm is violated by being too big
    idxs_J = (idxs_J_to_enf == idx_J) | (V > 1) ;
    
    % enforce the norm
    P_out(J,idxs_J) = P_out(J,idxs_J) ./ repmat(V(idxs_J),length(J),1) ;   
end