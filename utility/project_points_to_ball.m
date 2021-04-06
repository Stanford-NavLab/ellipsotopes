function P_out = project_points_to_ball(P_in,p_norm)
% P_out = project_points_to_ball(P_in,p_norm)
%
% Authors: Shreyas Kousik
% Created: 5 Apr 2021
% Updated: nop

    V = vecnorm(P_in,p_norm,1) ;
    P_out = P_in ./ repmat(V,size(P_in,1),1) ;
end