%% description
% This script tests computing the convex hull of two ellipsotopes.
%
% Authors: Adam Dai
% Created: 19 Apr 2021
% Updated: 
clear ; clc
%% user parameters

% ellipsotopes
p_norm = 2;
c1 = zeros(2,1);
G1 = 2*rand(2,4) - 1;
A1 = [1 0 -1 0; 0 1 0 -1];
b1 = [1; 0];

c2 = zeros(2,1);
G2 = 2*rand(2,4) - 1;
A2 = [1 0 -1 0; 0 1 0 -1];
b2 = [1; 0];

%% automated from here

% retrieve some dimensions
n = size(c1,1);
q1 = size(A1,1); m1 = size(A1,2);
q2 = size(A2,1); m2 = size(A2,2);

% create convex hull ellipsotope
c_CH = 0.5*(c1+c2);
G_CH = [G1, G2, 0.5*(c1-c2), zeros(n,2*(m1+m2))];

A_31 = [eye(m1); -eye(m1); zeros(m2,m1); zeros(m2,m1)];
A_32 = [zeros(m1,m2); zeros(m1,m2); eye(m2); -eye(m2)];
A_30 = [-0.5*ones(2*m1,1); 0.5*ones(2*m2,1)];

A_CH = [A1,           zeros(q1,m2), -b1/2, zeros(q1,2*(m1+m2));
        zeros(q2,m1), A2,            b2/2, zeros(q1,2*(m1+m2));
        A_31,         A_32,          A_30, eye(2*(m1+m2))];
    
b_CH = [b1/2; b2/2; -ones(2*(m1+m2),1)];
I_CH = {1:m1, m1+1:m1+m2, m1+m2+1:3*(m1+m2)+1};

%% plotting
figure(1); clf; axis equal; hold on; grid on;

% plot etopes
E1 = ellipsotope(p_norm,c1,G1,A1,b1);
E2 = ellipsotope(p_norm,c2,G2,A2,b2);
E_CH = ellipsotope(p_norm,c_CH,G_CH,A_CH,b_CH,I_CH);
plot(E1); plot(E2); plot(E_CH);


