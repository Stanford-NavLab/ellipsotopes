function heur = compute_order_reduction_MVOE_heuristic(G_i,G_j)
% heur = compute_order_reduction_MVOE_heuristic(G_i,G_j)
%
% Compute the heuristic value to estimate the volume of the Minkowski sum
% of a pair of ellipsoids represented as basic 2-ellipsotopes with
% generator matrices G_i and G_j
%
% Authors: Shreyas Kousik
% Created: 5 Mar 2022
% Updated: --
    Q_i = inv(pinv(G_i)'*pinv(G_i)) ;
    Q_j = inv(pinv(G_j)'*pinv(G_j)) ;
    Q_MVOE = 2*(Q_i + Q_j) ; % heuristic assumes bt = 1 ;
    heur = 1/det(inv(Q_MVOE)) ;
end