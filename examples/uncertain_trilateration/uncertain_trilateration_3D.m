%% description
% This script demonstrates an application of ellipsotopes for trilateration
% with uncertain range measurements and uncertain beacon positions. 
%
% Authors: Adam Dai
% Created: 22 Apr 2021
% Updated:

%% 3-D example

n = 3; % dimension
n_sats = 3;

% sphere of radius 20 representing satellite orbit, and sphere of radius 10
% representing Earth's surface
R = 20; r = 10;
[orbit_F,orbit_V] = make_sphere(R,zeros(3,1));
[earth_F,earth_V] = make_sphere(r,zeros(3,1));
figure(1) ; clf ; axis equal ; hold on ; grid on ; view(3) ;

% plot both spheres
patch('faces',orbit_F,'vertices',orbit_V,'facealpha',0.05,'edgealpha',0.2)
patch('faces',earth_F,'vertices',earth_V,'facealpha',0.1,'edgealpha',0.2,'facecolor','blue')

% satellite positions matrix
% (taken from orbit sphere)
S = zeros(n,n_sats); 
[S(1,1),S(2,1),S(3,1)] = sph2cart(0.5,0.4,R);
[S(1,2),S(2,2),S(3,2)] = sph2cart(0.1,-0.2,R);
[S(1,3),S(2,3),S(3,3)] = sph2cart(-0.4,0.2,R);
figure(1); scatter3(S(1,:),S(2,:),S(3,:),'bd','filled')

% satellite uncertainties (same for all satellites)
% represented with generator matrix
U_B = diag([0.1,0.2,0.1]);
 
% beacon uncertain position ellipsotopes
S1 = ellipsotope(2,S(:,1),U_B);
S2 = ellipsotope(2,S(:,2),U_B);
S3 = ellipsotope(2,S(:,3),U_B);

% true receiver position (taken from Earth surface)
x = zeros(n,1); r = 10;
[x(1),x(2),x(2)] = sph2cart(0.0,0.0,r);
figure(1); scatter3(x(1),x(2),x(3),'r*')

% (exact) range measurements
z = vecnorm(S - x)';

% measurement uncertainty (same for all measurements)
U_z_G = [0.2 0.1 0.3 
         0.1 0.2 0.1
         0.1 0.1 0.2];
U_z = ellipsotope(2,zeros(n,1),U_z_G);

% (exact) measurement ellipsotopes
Z1 = ellipsotope(2,zeros(n,1),z(1)*eye(n));
Z2 = ellipsotope(2,zeros(n,1),z(2)*eye(n));
Z3 = ellipsotope(2,zeros(n,1),z(3)*eye(n));

% + beacon uncertainty + measurement uncertainty
Z1 = Z1 + S1 + U_z;
Z2 = Z2 + S2 + U_z;
Z3 = Z3 + S3 + U_z;

% intersect uncertain measurement ellipsotopes
Z = Z1 & Z2 & Z3;

% plot
figure(1); axis equal; grid on; hold on
%plot(B1); plot(B2); plot(B3)
plot(Z1,'facealpha',0.1,'edgealpha',0.1); 
plot(Z2,'facealpha',0.1,'edgealpha',0.1); 
plot(Z3,'facealpha',0.1,'edgealpha',0.1);
plot(Z,'facecolor','r','edgecolor','r','facealpha',0.1);
