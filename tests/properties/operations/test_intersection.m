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
E1 = ellipsotope(2,[0;0],[1 0.5 -0.5; 0 0.866 0.866],[],[],{1,2,3});
E2 = ellipsotope(2,[0;0],diag([1;2]));

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
% perform the intersection 
E_int = E1 & E2;

%% plotting
h = figure(1) ; clf ; axis equal ; hold on ; grid on ;

plot(E1,'facecolor','r','edgecolor','r','facealpha',0.1);
plot(E2,'facecolor','b','edgecolor','b','facealpha',0.1);
plot(E_int,'facecolor','m','edgecolor','m','facealpha',1.0);
xlabel('$x\langle 1 \rangle$','Interpreter','latex'); ylabel('$x\langle 2 \rangle$','Interpreter','latex');

lim = axis; axis(lim + [-1 1 -1 1]);
legend('$E_1$','$E_2$','$E_1\cap E_2$','Interpreter','latex');

if flag_save_figure
    save_figure_to_pdf(h,'intersection.pdf')
    save_figure_to_png(h,'intersection.png')
end