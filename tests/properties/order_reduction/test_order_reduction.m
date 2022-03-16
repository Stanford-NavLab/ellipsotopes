%% description
% This script tests generator-based order reduction.
%
% Authors: Adam Dai and Shreyas Kousik
% Created: 20 May 2021
% Updated: 16 Mar 2022 (added rng seed)
%
%% user parameters
% rng seed
rng(1) ;

% basic 2-ellipsotope
% G = rand(2,5);
% E = ellipsotope(2,[0;0],G);

E = make_random_ellipsotope(2,2,20);
disp(['Original ellipsotope number of generators: ', num2str(E.order)]);

% whether or not to save flag
flag_save_figure = false ;

%% automated from here
% perform the order reduction
E_red = reduce(E,7);  % remove 7 generators
disp(['Reduced ellipsotope number of generators: ', num2str(E_red.order)]);

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