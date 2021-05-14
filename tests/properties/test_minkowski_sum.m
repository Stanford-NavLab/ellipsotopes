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
E2 = ellipsotope(2,[0;0],diag([1;2]));

% or generate random ellipsotopes
gen_random_flag = true;

% whether or not to save flag
flag_save_figure = true;

%% automated from here

if gen_random_flag
    E1 = make_random_ellipsotope();
    E2 = make_random_ellipsotope();
end

E_sum = E1 + E2;

%% plotting
fh = figure(1) ; clf ; axis equal ; hold on ; grid on ;
plot_ray_tracing_2D(E_sum,200,'facecolor','m','edgecolor','m','facealpha',0.1);
plot_ray_tracing_2D(E1,200,'facecolor','r','edgecolor','r','facealpha',0.1);
plot_ray_tracing_2D(E2,200,'facecolor','b','edgecolor','b','facealpha',0.1);
xlabel('$x\langle 1 \rangle$','Interpreter','latex'); ylabel('$x\langle 2 \rangle$','Interpreter','latex');
%set(gca,'fontsize',15)
lims = axis; axis(1.2*lims);

if flag_save_figure
    save_figure_to_pdf(fh,'minkowski_sum.pdf')
end