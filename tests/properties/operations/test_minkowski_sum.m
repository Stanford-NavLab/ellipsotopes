%% description
% This script tests computing the Minkowski sum of two ellipsotopes.
%
% Authors: Adam Dai 
% Created: 15 Mar 2021
% Updated: 13 May 2021
%
%% user parameters
% specify the two ellipsotopes
E1 = ellipsotope(2,[0;0],[1 0.5 -0.5; 0 0.866 0.866],[],[],{1,2,3});
E2 = rotation_matrix_2D(0.2) * ellipsotope(2,[0;0],diag([2;1]));

% or generate random ellipsotopes
gen_random_flag = false;

% whether or not to save flag
flag_save_figure = true;

%% automated from here
% generate random input if desired
if gen_random_flag
    E1 = make_random_ellipsotope();
    E2 = make_random_ellipsotope();
end
% perform the minkowski sum
E_sum = E1 + E2;

%% plotting
c1 = [0, 0.4470, 0.7410];
c2 = [0.8500, 0.3250, 0.0980];
c3 = [0.4940, 0.1840, 0.5560];
h = figure(1) ; clf ; axis equal ; hold on ; grid on ;
plot(E1,'facecolor',c1,'edgecolor',c1,'facealpha',0.1);
plot(E2,'facecolor',c2,'edgecolor',c2,'facealpha',0.1);
plot(E_sum,'facecolor',c3,'edgecolor',c3,'facealpha',0.1);
xlabel('$x\langle 1 \rangle$','Interpreter','latex'); ylabel('$x\langle 2 \rangle$','Interpreter','latex');
%set(gca,'fontsize',15)
lim = axis; axis(lim + [-1 1 -1 1]);
legend('$E_1$','$E_2$','$E_1\oplus E_2$','Interpreter','latex');

if flag_save_figure
    save_figure_to_pdf(h,'minkowski_sum.pdf')
    save_figure_to_png(h,'minkowski_sum.png')
end