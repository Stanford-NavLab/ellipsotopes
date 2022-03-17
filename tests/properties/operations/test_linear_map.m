%% description
% Test linear map of an ellipsotope.
%
% Authors: Adam Dai 
% Created: 19 May 2021
% Updated: 
%
%% user parameters
% specify the ellipsotopes and map
E = ellipsotope(2,[0;0],[1 0.5; 0 0.866],[],[],{1,2});
%E = E + ellipsotope(2,[0;0],eye(2));
A = [1 0.8; 0.4 2];

% or generate random ellipsotopes
gen_random_flag = false;

% whether or not to save flag
flag_save_figure = true;

%% automated from here
% generate random input if desired
if gen_random_flag
    E = make_random_ellipsotope();
end
% perform the linear map 
E_lm = A * E;

%% plotting
h = figure(1) ; clf ; axis equal ; hold on ; grid on ;

plot(E,'facecolor','r','edgecolor','r','facealpha',0.1);
plot(E_lm,'facecolor','b','edgecolor','b','facealpha',0.1);
xlabel('$x\langle 1 \rangle$','Interpreter','latex'); ylabel('$x\langle 2 \rangle$','Interpreter','latex');

lim = axis; axis(lim + [-1 1 -1 1]);
legend('$E$','$AE$','Interpreter','latex');

if flag_save_figure
    save_figure_to_pdf(h,'linear_map.pdf')
    save_figure_to_png(h,'linear_map.png')
end