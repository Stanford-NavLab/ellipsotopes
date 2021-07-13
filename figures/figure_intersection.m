%% description
% This script intersects 2 circular ellipsotopes to produce a constrained
% ellipsotope
%
% Authors: Adam Dai and Shreyas Kousik
% Created: shrug
% Updated: 9 July 2021 (resized stuff for prettiness)
%
%% user parameters
% 2 ellipses
E1 = ellipsotope(2,[0;0],diag([2;1]));
E2 = ellipsotope(2,[1.5;0],diag([1;2]));

% general ellipsotopes
rng(5)
E3 = make_random_ellipsotope(2,2,5,0,3) ;
E4 = make_random_ellipsotope(2,2,3,0,3) ;

%% automated from here
% perform the intersection 
E_int1 = E1 & E2;
E_int2 = E3 & E4;

%% plotting setup
% colors
c1 = [0.7 0 0.7];
c2 = [0 0.7 0.7];
c3 = [0.5 0.5 0];

%% plotting
% h_1 = figure(1) ; clf ; axis equal ; hold on ; grid on ;
% plot(E1,'color',c1,'facealpha',0.1,'edgealpha',0.7);
% plot(E2,'color',c2,'facealpha',0.1,'edgealpha',0.7);
% plot(E_int1,'color',c3,'facealpha',0.1,'edgealpha',1);
% % xlabel('$x\langle 1 \rangle$','Interpreter','latex'); ylabel('$x\langle 2 \rangle$','Interpreter','latex');
% % lim = axis; axis(lim + [-1 1 -1 1]);
% axis tight
% legend('$E_1$','$E_2$','$E_1\cap E_2$','Interpreter','latex');
% 
% % save_figure_to_pdf(h,'intersection_ellipses.pdf')
% save_figure_to_png(h_1,'intersection_ellipses.png')

%% figure 2
h_2 = figure(2) ; clf ; axis equal ; hold on ; grid on ;
plot(E3,'color',c1,'facealpha',0.1,'edgealpha',0.7);
plot(E4,'color',c2,'facealpha',0.1,'edgealpha',0.7);
plot(E_int2,'color',c3,'facealpha',0.1,'edgealpha',1);
% xlabel('$x\langle 1 \rangle$','Interpreter','latex'); ylabel('$x\langle 2 \rangle$','Interpreter','latex');
% lim = axis; axis(lim + [-1 1 -1 1]);
legend('$E_1$','$E_2$','$E_1\cap E_2$','Interpreter','latex');
set_plot_fontsize(15)

% save_figure_to_pdf(h,'intersection_general.pdf')
save_figure_to_png(h_2,'intersection_general.png')