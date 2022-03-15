%% description
% This script tests overapproximation of a body rotated by an uncertain
% amount using zonotopes
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
Z_body = zonotope(c,G);

h_sigma = 0.1;
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
Z_rot1 = rotation_matrix_2D(d_h) * Z_body;
Z_rot2 = rotation_matrix_2D(-d_h) * Z_body;

% outer circle
E_circ = ellipsotope(2,zeros(2,1),r*eye(2));

% bounding dimensions
z = (w/2)*sin(pi/2-d_h) + (l/2)*sin(d_h);
psi = atan2(w,l);
if d_h < psi
    H = r * cos(psi - d_h);
else
    H = r;
end

% creating bounding zonotope
G_bound = diag([H,z]);
Z_bound = rotation_matrix_2D(theta) * zonotope(c,G_bound);

%% plot
figure(1); hold on; grid on; axis equal
plot(Z_bound,[1,2],'Filled',true,'FaceColor','r','EdgeColor','r','FaceAlpha',0.5,'EdgeAlpha',1.0);
plot(Z_rot1,[1,2],'Filled',true,'FaceAlpha',0.1,'FaceColor','b','EdgeColor','b','LineWidth',1.5); 
plot(Z_rot2,[1,2],'Filled',true,'FaceAlpha',0.1,'FaceColor','b','EdgeColor','b','LineWidth',1.5); 

lim = axis; axis(lim + 0.25*[-1 1 -1 1]);
legend('Overbounding of \newline heading uncertainty','Rotated body');
set(gca,'fontsize',15)