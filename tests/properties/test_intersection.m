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
% E1 = ellipsotope(2,[0;0],diag([1,2]));
% E2 = ellipsotope(2,[1;0],eye(2));

p = 2 ;
c = zeros(2,1) ;
G = 2*rand(2,4) - 1;
A = [-1 1 -1 1] ;
b = 0.5 ;
I = {[1,2],[3,4]} ;
E1 = ellipsotope(p,c,G,A,b,I);

c = ones(2,1);
G = 2*rand(2,4) - 1;
E2 = ellipsotope(p,c,G,A,b,I);

%% automated from here
% perform the intersection 
E_int = E1 & E2;

%% plotting
h = figure(1) ; clf ; axis equal ; hold on ; grid on ;

plot(E1);
plot(E2);
plot(E_int,'facecolor','r','edgecolor','r','facealpha',0.1);

set(gca,'fontsize',15)