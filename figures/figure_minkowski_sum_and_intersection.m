%% description
% This script generates a figure illustrating minkowski sum and 
% intersection of ellipsotopes.
%
% Authors: Adam Dai 
% Created: 15 Mar 2021
% Updated: 9 July 2021 (Shreyas fiddled with more things)
%
%% minkowski sum
% rng seed
rng(0)

% specify the two ellipsotopes
% E1 = ellipsotope(2,[0;0],[1 0.5 -0.5; 0 0.866 0.866],[],[],{1,2,3});
E1 = make_random_ellipsotope(2,2,3,0,3) ;
E1 = E1 + [-1; -1];
E1 = rotation_matrix_2D(0.4) * E1;
E2 = rotation_matrix_2D(0.2) * ellipsotope(2,[0;0],diag([2;1]));

% perform the minkowski sum
E_sum = E1 + E2;

% perform the intersection 
E_int = E1 & E2;

% %% intersection
% % general ellipsotopes
% rng(5)
% E3 = make_random_ellipsotope(2,2,5,0,3) ;
% E4 = make_random_ellipsotope(2,2,3,0,3) ;
% 
% % perform the intersection 
% E_int = E3 & E4;

%% plotting
% colors
c1 = [0.7 0 0.7];
c2 = [0 0.2 0.7];
c3 = [0.5 0.5 0];
c4 = [0.7 0 0];

h = figure(1) ; clf ; axis equal ; hold on ; grid on ; axis tight
plot(E1,'color',c1,'facealpha',0.1);
plot(E2,'color',c2,'facealpha',0.1);
plot(E_sum,'color',c3,'facealpha',0.1);
plot(E_int,'color',c4,'facealpha',0.1,'edgealpha',1,'num_points',500);
legend('$E_1$','$E_2$','$E_1\oplus E_2$','$E_1\cap E_2$','Interpreter','latex');

% % plot intersection
% subplot(1,2,2);
% plot(E3,'color',c1,'facealpha',0.1,'edgealpha',0.7);
% plot(E4,'color',c2,'facealpha',0.1,'edgealpha',0.7);
% plot(E_int2,'color',c3,'facealpha',0.1,'edgealpha',1);
% legend('$E_1$','$E_2$','$E_1\cap E_2$','Interpreter','latex');


%% labeling
% xlabel('$x\langle 1 \rangle$','Interpreter','latex'); ylabel('$x\langle 2 \rangle$','Interpreter','latex');
% lim = axis; axis(lim + [-1 1 -1 1]);
lim = axis; axis(lim + 0.25*[-2 2 -1 1]);
%set_plot_fontsize(15)

save_figure_to_pdf(h,'minkowski_sum_and_intersection.pdf')
save_figure_to_png(h,'minkowski_sum_and_intersection.png')