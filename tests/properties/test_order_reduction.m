%% description
% This script tests computing the Minkowski sum of two ellipsotopes.
%
% Authors: Adam Dai 
% Created: 20 May 2021
% Updated: 
%
%% user parameters
% basic 2-ellipsotope
% G = rand(2,5);
% E = ellipsotope(2,[0;0],G);

E = make_random_ellipsotope(2,2,20);

% whether or not to save flag
flag_save_figure = true;

%% automated from here
% perform the order reduction
E_red = reduce(E,13);

%% plotting
h = figure(1) ; clf ; axis equal ; hold on ; grid on ;
plot(E,'facecolor','r','edgecolor','r','facealpha',0.1);
plot(E_red,'facecolor','b','edgecolor','b','facealpha',0.1);
xlabel('$x\langle 1 \rangle$','Interpreter','latex'); ylabel('$x\langle 2 \rangle$','Interpreter','latex');
%set(gca,'fontsize',15)
lim = axis; axis(lim + [-1 1 -1 1]);
legend('$E$','$E$ reduced','Interpreter','latex');

if flag_save_figure
    save_figure_to_pdf(h,'order_reduction.pdf')
    save_figure_to_png(h,'order_reduction.png')
end