%% description
% This script tests plotting ellipsotopes with more than 3 generators
%
% Authors: Shreyas Kousik
% Created: 15 Mar 2021
%
%% user parameters
% rng seed
rng(10) ;

% p norm
p = 10 ;

% center
c = rand(2,1) ;

% generators
G = 2*rand(2,2) - 1 ;

%% automated from here
% make ellipsotope
E = ellipsotope(p,c,G) ;

% make comparable zonotope (requires CORA)
Z = zonotope(c,G) ;

%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ;

plot(Z) ;
plot(E) ;