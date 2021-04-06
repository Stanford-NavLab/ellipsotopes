function P_out = project_points_to_ball_product(P_in,p_norm,I,idxs_J_to_enf)
% P_out = project_points_to_ball_product(P_in,p_norm,I)
% P_out = project_points_to_ball_product(P_in,p_norm,I,idxs_J_to_enf)
%
% See also: test_projection_to_generalized_etope_ball
%
% Authors: Shreyas Kousik
% Created: 5 Apr 2021
% Updated: nup

% set default args in
if nargin < 4
    % make random list of which norm to enforce per point
    idxs_J_to_enf = rand_int(1,n_I,[],[],1,n_P) ;
end

% number of index subsets
n_I = length(I) ;

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