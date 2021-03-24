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

% initial e-tope
c = [0; 0]; % center
G = [1 0 2;
     0 1 -1]; % generator matrix 
 
% system
A = [1 0.1;
     0.1 1];
B = [0.5 0;
     0 0.5];

% nominal trajectory
N = 10;
U = ones(2,N);

% norm to consider
p = 2 ;

x_0 = ellipsotope(p,c,G);

%% automated from here
% reachability
R = {};

x = x_0;
for i = 1:N
    u_nom = U(:,i);
    u = u_nom;
    x = A * x + B * u; 
    R{i} = x;
end

%% plotting
for i = 1:N
   plot(R{i}); 
end