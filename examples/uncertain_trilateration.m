%% description
% This script demonstrates an application of ellipsotopes for trilateration
% with uncertain range measurements and uncertain beacon positions. 
%
% Authors: Adam Dai
% Created: 22 Apr 2021
% Updated:

%% 2-D example

% beacon positions matrix
B = [0 10 5;
     0 0  10];

% beacon uncertainties (same for all beacons)
% represented with generator matrix
U_B = diag([0.1,0.2]);
 
% beacon uncertain position ellipsotopes
B1 = ellipsotope(2,B(:,1),U_B);
B2 = ellipsotope(2,B(:,2),U_B);
B3 = ellipsotope(2,B(:,3),U_B);

% true receiver position
x = [7;3];

% (exact) range measurements
z = vecnorm(B - x)';

% measurement uncertainty (same for all measurements)
U_z_G = [0.2 0.1 
         0.1 0.2];
U_z = ellipsotope(2,zeros(2,1),U_z_G,[],[],{1,2});

% (exact) measurement ellipsotopes
Z1 = ellipsotope(2,zeros(2,1),z(1)*eye(2));
Z2 = ellipsotope(2,zeros(2,1),z(2)*eye(2));
Z3 = ellipsotope(2,zeros(2,1),z(3)*eye(2));

% + beacon uncertainty + measurement uncertainty
Z1 = Z1 + B1 + U_z;
Z2 = Z2 + B2 + U_z;
Z3 = Z3 + B3 + U_z;

% intersect uncertain measurement ellipsotopes
Z = Z1 & Z2 & Z3;

% plot
figure(); axis equal; grid on
plot(B1); plot(B2); plot(B3)
plot(Z1); plot(Z2); plot(Z3)
plot(Z,'facecolor','r','edgecolor','r','facealpha',0.1);
