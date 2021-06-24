%% description
% This script tests computing the Minkowski sum of two ellipsotopes.
%
% Authors: Adam Dai 
% Created: 15 Mar 2021
% Updated: 15 Jun 2021 (Shreyas fiddled with things)
%
%% user parameters
% rng seed
rng(0)

% specify the two ellipsotopes
% E1 = ellipsotope(2,[0;0],[1 0.5 -0.5; 0 0.866 0.866],[],[],{1,2,3});
E1 = make_random_ellipsotope(2,2,3,0,3) ;
E2 = rotation_matrix_2D(0.2) * ellipsotope(2,[0;0],diag([2;1]));

%% automated from here
% perform the minkowski sum
E_sum = E1 + E2;

%% plotting
c1 = [0.7 0 0.7];
c2 = [0 0.7 0.7];
c3 = [0.5 0.5 0];

h = figure(1) ; clf ; axis equal ; hold on ; grid on ;
plot(E1,'color',c1,'facealpha',0.1);
plot(E2,'color',c2,'facealpha',0.1);
plot(E_sum,'color',c3,'facealpha',0.1);
xlabel('$x\langle 1 \rangle$','Interpreter','latex'); ylabel('$x\langle 2 \rangle$','Interpreter','latex');
%set(gca,'fontsize',15)
lim = axis; axis(lim + [-1 1 -1 1]);
legend('$E_1$','$E_2$','$E_1\oplus E_2$','Interpreter','latex');

save_figure_to_pdf(h,'minkowski_sum.pdf')
save_figure_to_png(h,'minkowski_sum.png')