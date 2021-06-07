%% description
% This script tests a simple reachability example using ellipsotopes.
% An initial ellipsotope is created and propagated through a linear 
% dynamical system.
%
% Authors: Adam Dai
% Created: 24 May 2021
% Updated: 4 June 2021
%
clc
%% user parameters

theta = 2*pi/3;
w = 1; l = 2;
r = norm([w l]) / 2;
c = zeros(2,1);
G = rotation_matrix_2D(theta) * 0.5 * diag([l w]);
E_body = ellipsotope(2,zeros(2,1),G,[],[],{1,2});

h_sigma = 0.01;
P = 0.9;
d_h = erfinv(P) * h_sigma * sqrt(2);

% sanity check heading confidence interval
N_samples = 1000;
h_samples = normrnd(0,h_sigma,N_samples,1);
in_count = 0;
for i = 1:N_samples
    if abs(h_samples(i)) < d_h
        in_count = in_count + 1;
    end
end
disp(['P(inside):', num2str(in_count/N_samples)]);

% rotated body
E_rot1 = rotation_matrix_2D(d_h) * E_body;
E_rot2 = rotation_matrix_2D(-d_h) * E_body;

% outer circle
E_circ = ellipsotope(2,zeros(2,1),r*eye(2));

% halfplanes
z = (w/2)*sin(pi/2-d_h) + (l/2)*sin(d_h);
x1 = c + z*[cos(pi/2+theta);sin(pi/2+theta)];
h1 = [cos(theta+pi/2),sin(theta+pi/2)]; f1 = h1*x1;
x2 = c - z*[cos(pi/2+theta);sin(pi/2+theta)];
h2 = -[cos(theta+pi/2),sin(theta+pi/2)]; f2 = h2*x2;

% intersect
E_rot = halfspace_intersect(E_circ,h1,f1);
E_rot = halfspace_intersect(E_rot,h2,f2);

%% plot
figure(1); hold on; grid on; axis equal
plot(E_circ,'EdgeAlpha',1.0,'FaceAlpha',0.0,'EdgeColor','#D95319','LineWidth',1.0,'LineStyle','--');
plot(E_rot,'FaceColor','r','EdgeColor','r','FaceAlpha',0.5,'EdgeAlpha',1.0);
plot(E_rot1,'FaceAlpha',0.1,'FaceColor','b','EdgeColor','b','LineWidth',1.5); 
plot(E_rot2,'FaceAlpha',0.1,'FaceColor','b','EdgeColor','b','LineWidth',1.5); 

lim = axis; axis(lim + 0.25*[-1 1 -1 1]);
legend('Circumscribing circle','Overbounding of \newline heading uncertainty','Rotated body');
set(gca,'fontsize',15)