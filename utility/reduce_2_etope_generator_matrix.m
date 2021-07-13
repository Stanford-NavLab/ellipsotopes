function [G_out, Q] = reduce_2_etope_generator_matrix(G)
% G = reduce_2_etope_generator_matrix(G)
% [G,Q] = reduce_2_etope_generator_matrix(G)
%
% Reduce the generator matrix of a basic ellipsotope. Can also return the
% ellipsoid shape matrix Q >= 0 as a second output.
%
% Note: Works for "wide" (fat) generator matrices.
%
% Authors: Shreyas Kousik
% Created: 10 Feb 2021
% Updated: 13 Jul 2021
G_Q = pinv(G) ;
Q = G_Q'*G_Q ;
G_out = pinv(Q^(1/2)) ;
end