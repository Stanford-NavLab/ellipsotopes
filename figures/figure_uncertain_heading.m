%% description
% This script plots the heading uncertainty ellipsotope along with the
% corresponding extreme rotated bodies.
%
% Authors: Adam Dai
% Created: 24 May 2021
% Updated: 17 Mar 2022 (extended the halfspace lines, added plot flag)
%
clear ; clc ; close all
%% user parameters

theta = pi/6; %2*pi/3 ;
w = 2; l = 4;
r = norm([w l]) / 2;
c = zeros(2,1);
G = rotation_matrix_2D(theta) * 0.5 * diag([l w]);
E_body = ellipsotope(2,zeros(2,1),G,[],[],{1,2});

h_sigma = 0.1;
P = 0.9;
d_h = erfinv(P) * h_sigma * sqrt(2);

% sanity check heading confidence interval
N_samples = 1000;

% whether or not to save figure
flag_save_figure = true ;

%% automated from here
% test how many points are inside the tope, given the confidence bound on
% the Gaussian distribution
h_samples = normrnd(0,h_sigma,N_samples,1);
in_count = 0;
for i = 1:N_samples
    if abs(h_samples(i)) < d_h
        in_count = in_count + 1;
    end
end
disp(['P(inside): ', num2str(in_count/N_samples)]);

%% plotting setup
% rotated body
E_rot1 = rotation_matrix_2D(d_h) * E_body;
E_rot2 = rotation_matrix_2D(-d_h) * E_body;

% outer circle
E_circ = ellipsotope(2,zeros(2,1),r*eye(2));

% halfspaces
z = (w/2)*sin(pi/2-d_h) + (l/2)*sin(d_h);
x1 = c + z*[cos(pi/2+theta);sin(pi/2+theta)];
h1 = [cos(theta+pi/2),sin(theta+pi/2)]; f1 = h1*x1;
x2 = c - z*[cos(pi/2+theta);sin(pi/2+theta)];
h2 = -[cos(theta+pi/2),sin(theta+pi/2)]; f2 = h2*x2;

% intersect
E_rot = halfspace_intersect(E_circ,h1,f1);
E_rot = halfspace_intersect(E_rot,h2,f2);

% clip top and bottom
psi = atan2(w,l);
if d_h < psi
    H = r * cos(psi - d_h);
    x3 = c + H * [cos(theta);sin(theta)];
    h3 = [cos(theta),sin(theta)]; 
    f3 = h3 * x3;
    x4 = c - H * [cos(theta);sin(theta)];
    h4 = -[cos(theta),sin(theta)]; 
    f4 = h4 * x4;
    E_rot = halfspace_intersect(E_rot,h3,f3);
    E_rot = halfspace_intersect(E_rot,h4,f4);
end

%% plotting
fh = figure(1); clf ; hold on ; grid on ; axis equal

plot(E_circ,'EdgeAlpha',1.0,'FaceAlpha',0.0,'EdgeColor','#D95319','LineWidth',1.0);%,'LineStyle','--');
plot(E_rot,'FaceColor',[0.5 0.5 0],'EdgeColor',[0.5 0.5 0],'FaceAlpha',0.5,'EdgeAlpha',1.0,'num_points',10000);
%plot(E_rot,'FaceColor',[0.5 0.5 0],'EdgeColor',[0.5 0.5 0],'FaceAlpha',0.5,'EdgeAlpha',1.0,'num_points',500);
plot(E_rot1,'FaceAlpha',1.0,'FaceColor',[0.7 0.7 1],'EdgeColor','b','LineWidth',1.5); 
plot(E_body,'FaceAlpha',1.0,'FaceColor',[0.7 0.7 1],'EdgeColor','b','LineWidth',1.5); 
plot(E_rot2,'FaceAlpha',1.0,'FaceColor',[0.7 0.7 1],'EdgeColor','b','LineWidth',1.5); 

% show halfspaces
plot([x1(1) - 4*h1(2), x1(1) + 4*h1(2)],[x1(2) + 4*h1(1), x1(2) - 4*h1(1)],'--','color',[0 0.8 0.6],'LineWidth',2)
plot([x2(1) - 4*h2(2), x2(1) + 4*h2(2)],[x2(2) + 4*h2(1), x2(2) - 4*h2(1)],'--','color',[0 0.8 0.6],'LineWidth',2)
plot([x3(1) - 4*h3(2), x3(1) + 4*h3(2)],[x3(2) + 4*h3(1), x3(2) - 4*h3(1)],'--','color',[0 0.8 0.6],'LineWidth',2)
plot([x4(1) - 4*h4(2), x4(1) + 4*h4(2)],[x4(2) + 4*h4(1), x4(2) - 4*h4(1)],'--','color',[0 0.8 0.6],'LineWidth',2)

% marker for heading
marker_verts = 0.2 * rotation_matrix_2D(theta-d_h) * [2 -0.5 -0.5; 0 sqrt(3)/2 -sqrt(3)/2];
patch(marker_verts(1,:),marker_verts(2,:),'black');

% lim = axis; axis(lim + 0.25*[-2 2 -1 1]);
axis([-4,4,-2.4,2.4])
legend('Circumscribing Circle','Overbound Ellipsotope','Rotated Body',...
    'location','northwest');
set(gca,'fontsize',12)

if flag_save_figure
    % exportgraphics(fh,'heading_uncertainty_cropped.png') ;
    exportgraphics(fh,'heading_uncertainty_cropped.pdf') ;
end