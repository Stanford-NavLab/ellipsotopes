%% description
% This script plots several example ellipsotopes
%
% Authors: Shreyas Kousik
% Created: 15 Mar 2021
% Updated: 22 Mar 2021
%
%% zonotope
c = [0;0];
G = eye(2);
E_zono = ellipsotope(2,c,G,[],[],{1,2});

h = figure(1); clf; axis equal; hold on; grid on;
plot(E_zono);
xlabel('$x\langle 1 \rangle$','Interpreter','latex'); 
ylabel('$x\langle 2 \rangle$','Interpreter','latex');
lim = axis; axis(lim + 0.25*[-1 1 -1 1]);

% show generators
n_gen = size(G,2);
quiver(repmat(c(1),n_gen,1),repmat(c(2),n_gen,1),G(1,:)',G(2,:)','off')

%% ellipsoid
c = [0;0];
G = eye(2);
E_ellip = ellipsotope(2,c,G);

h = figure(1); clf; axis equal; hold on; grid on;
plot(E_ellip);
xlabel('$x\langle 1 \rangle$','Interpreter','latex'); 
ylabel('$x\langle 2 \rangle$','Interpreter','latex');
lim = axis; axis(lim + 0.25*[-1 1 -1 1]);

% show generators
n_gen = size(G,2);
quiver(repmat(c(1),n_gen,1),repmat(c(2),n_gen,1),G(1,:)',G(2,:)','off')
