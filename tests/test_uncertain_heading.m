%% description
% This script tests a simple reachability example using ellipsotopes.
% An initial ellipsotope is created and propagated through a linear 
% dynamical system.
%
% Authors: Adam Dai
% Created: 24 Mar 2021
% Updated: 
%
%% user parameters

w = 1; l = 2;
r = norm([w l]) / 2;
G = 0.5 * diag([w l]);
E_body = ellipsotope(2,zeros(2,1),G,[],[],{1,2});

% h_sigma = 0.1;
% P = 0.9;
% eps = -log((1-P)^2 * (h_sigma/(2*pi)));
% d_h = sqrt(eps*h_sigma);
d_h = 0.01;

% sanity check heading confidence interval
% N_samples = 1000;
% h_samples = normrnd(0,h_sigma,N_samples,1);
% in_count = 0;
% for i = 1:N_samples
%     if abs(h_samples(i)) < d_h
%         in_count = in_count + 1;
%     end
% end
% disp(['P(inside):', num2str(in_count/N_samples)]);

% rotated body
E_rot1 = rotation_matrix_2D(d_h) * E_body;
E_rot2 = rotation_matrix_2D(-d_h) * E_body;

% circle
E_circ = ellipsotope(2,zeros(2,1),r*eye(2));



%% plot
figure(1); hold on; grid on; axis equal
plot(E_body);
plot(E_rot1,'edgecolor','r','facecolor','r'); plot(E_rot2,'edgecolor','r','facecolor','r');
plot(E_circ,'edgecolor','m','facecolor','m');
