%% description
% This script tests computing the convex hull of two ellipsotopes.
%
% Authors: Adam Dai (and Shreyas sneakin in here too)
% Created: 19 Apr 2021
% Updated: 20 Apr 2021
clear ; clc
%% user parameters

% ellipsotopes
p_norm = 2;
c_1 = zeros(2,1);
G_1 = 2*rand(2,4) - 1;
A_1 = [1 0 -1 0; 0 1 0 -1];
b_1 = [1; 0];

c_2 = zeros(2,1);
G_2 = 2*rand(2,4) - 1;
A_2 = [1 0 -1 0; 0 1 0 -1];
b_2 = [1; 0];

% c_1 = zeros(2,1) ;
% G_1 = eye(2) ;
% A_1 = [] ;
% b_1 = [] ;
% 
% c_2 = [1;0] ;
% G_2 = eye(2) ;
% A_2 = [] ;
% b_2 = [] ;

%% automated from here
% retrieve some dimensions
n = size(c_1,1);
q_1 = size(A_1,1); m_1 = size(A_1,2);
q_2 = size(A_2,1); m_2 = size(A_2,2);

% create convex hull ellipsotope
c_CH = 0.5*(c_1+c_2);
G_CH = [G_1, G_2, 0.5*(c_1-c_2), zeros(n,2*(m_1+m_2))];

A_31 = [eye(m_1); -eye(m_1); zeros(m_2,m_1); zeros(m_2,m_1)];
A_32 = [zeros(m_1,m_2); zeros(m_1,m_2); eye(m_2); -eye(m_2)];
A_30 = [-0.5*ones(2*m_1,1); 0.5*ones(2*m_2,1)];

A_CH = [A_1,           zeros(q_1,m_2), -b_1/2, zeros(q_1,2*(m_1+m_2));
        zeros(q_2,m_1), A_2,            b_2/2, zeros(q_1,2*(m_1+m_2));
        A_31,         A_32,          A_30, eye(2*(m_1+m_2))];
    
b_CH = [b_1/2; b_2/2; -ones(2*(m_1+m_2),1)];
I_CH = {1:m_1, m_1+1:m_1+m_2, m_1+m_2+1:3*(m_1+m_2)+1};

%% plotting
figure(1); clf; axis equal; hold on; grid on;

% plot etopes
E1 = ellipsotope(p_norm,c_1,G_1,A_1,b_1);
E2 = ellipsotope(p_norm,c_2,G_2,A_2,b_2);
E_CH = ellipsotope(p_norm,c_CH,G_CH,A_CH,b_CH,I_CH);
plot(E1); plot(E2); plot(E_CH);


