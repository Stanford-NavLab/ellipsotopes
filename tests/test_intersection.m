%% description
% This script intersects 2 circular ellipsotopes to produce a constrained
% ellipsotope
%
% Authors: Adam Dai and Shreyas Kousik
% Created: shrug
% Updated: 15 Mar 2021
%
%% user parameters
% specify the two ellipsotopes
E1 = ellipsotope(2,[0;0],eye(2));
E2 = ellipsotope(2,[1;0],eye(2));

%% automated from here
% perform the intersection 
E_int = E1 & E2;

%% plotting
h = figure(1) ; clf ; axis equal ; hold on ; grid on ;

plot(E1);
plot(E2);
plot(E_int,'facecolor','r','edgecolor','r','facealpha',0.1);

set(gca,'fontsize',15)