%% description
% This script tests computing the convex hull of two ellipsotopes.
%
% Authors: Adam Dai (and Shreyas sneakin in here too)
% Created: 19 Apr 2021
% Updated: 15 June 2021 (changed plot slightly)
clear ; clc
%% user parameters
% rng seed
rng(0)

% specify ellipsotopes
p_norm = 2 ;

c_1 = zeros(2,1) ;
G_1 = eye(2) ;
A_1 = [] ;
b_1 = [] ;

c_2 = 2*rand(2,1) - 1 ;
G_2 = 2*rand(2) - 1 ;
A_2 = [] ;
b_2 = [] ;

%% automated from here
E1 = ellipsotope(p_norm,c_1,G_1,A_1,b_1);
E2 = ellipsotope(p_norm,c_2,G_2,A_2,b_2);

E_CH = convhull(E1,E2);

%% plotting
h = figure(1); clf; axis equal; hold on; grid on;

c1 = [0.7 0 0.7];
c2 = [0 0.7 0.7];
c3 = [0.5 0.5 0];

% plot etopes
plot(E1,'color',c1) ;
plot(E2,'color',c2) ;
plot(E_CH,'color',c3) ;
xlabel('$x\langle 1 \rangle$','Interpreter','latex'); ylabel('$x\langle 2 \rangle$','Interpreter','latex');
lim = axis; axis(lim + 0.25*[-1 1 -1 1]);
legend('$E_1$','$E_2$','$E_{CH}$','Interpreter','latex');

save_figure_to_pdf(h,'convex_hull.pdf')
save_figure_to_png(h,'convex_hull.png')