%% description
% This script is for plotting example ellipsotopes
%
% Authors: Adam Dai
% Created: 24 May 2021
% Updated: 
%

% robot body 
G_robot = rotation_matrix_2D(pi*0.2) * diag([1 2]);
E_robot = ellipsotope(2,zeros(2,1),G_robot,[],[],{1,2});

% robot uncertainty
G_uncrt = rotation_matrix_2D(pi*-0.2) * diag([1 2]);
E_uncrt = ellipsotope(2,zeros(2,1),G_uncrt);

% mink sum
E_mink = E_robot + E_uncrt;

f = figure(2); 
subplot(1,2,1); axis equal; grid on;
plot(E_uncrt,'facecolor','r','edgecolor','r','facealpha',0.4)
xlabel('Uncertainty in $p\langle 1 \rangle$','Interpreter','latex'); ylabel('Uncertainty in $p\langle 2 \rangle$','Interpreter','latex');
lim = axis; axis(lim + [-1 1 -1 1]);

subplot(1,2,2); axis equal; grid on;
plot(E_robot,'facecolor','b','edgecolor','b','facealpha',1.0)
plot(E_mink,'facecolor','m','edgecolor','m','facealpha',0.4)
xlabel('$p\langle 1 \rangle$','Interpreter','latex'); ylabel('$p\langle 2 \rangle$','Interpreter','latex');
legend('Robot body','Sum of body and uncertainty');
lim = axis; axis(lim + [-1 1 -1 1]);

save_figure_to_pdf(f,'path_planning_uncertainty_sum.pdf')
