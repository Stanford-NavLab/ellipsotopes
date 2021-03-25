%% Description
% Minkowski sum of 2 circles

%% user parameters
% define the ellipsotopes to be summed

c1 = [0;0]; G1 = [1 0 0 0; 0 1 0 0];
J1 = {[1,2],[3,4]}; A1 = [1 0 -1 0; 0 1 0 -1]; b1 = [1; 0];
E1 = ellipsotope(2, c1, G1, A1, b1, J1);

c2 = [1;1]; G2 = [1 0 0 0; 0 2 0 0];
J2 = {[1,2],[3,4]}; A2 = [1 0 -1 0; 0 2 0 -1]; b2 = [1; 0];
E2 = ellipsotope(2, c2, G2, A2, b2, J2);

E_sum = E1 + E2;


figure(1) ; clf ; axis equal ; hold on ; grid on ;

plot(E_sum,'facecolor','m','edgecolor','m','facealpha',0.1);
plot(E1,'facecolor','r','edgecolor','r','facealpha',0.1);
plot(E2,'facecolor','b','edgecolor','b','facealpha',0.1);
