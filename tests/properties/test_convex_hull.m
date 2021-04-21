%% description
% This script tests computing the convex hull of two ellipsotopes.
%
% Authors: Adam Dai (and Shreyas sneakin in here too)
% Created: 19 Apr 2021
% Updated: 20 Apr 2021
clear ; clc
%% user parameters
% rng seed
rng(0)

% ellipsotopes
p_norm = 2;
% c_1 = zeros(2,1);
% G_1 = 2*rand(2,4) - 1;
% A_1 = [1 0 -1 0; 0 1 0 -1];
% b_1 = [1; 0];
% 
% c_2 = zeros(2,1);
% G_2 = 2*rand(2,4) - 1;
% A_2 = [1 0 -1 0; 0 1 0 -1];
% b_2 = [1; 0];

c_1 = zeros(2,1) ;
G_1 = eye(2) ;
A_1 = [] ;
b_1 = [] ;

c_2 = [1;1] ;
G_2 = eye(2) ;
A_2 = [] ;
b_2 = [] ;

%% automated from here
% retrieve sizes of things
n = size(c_1,1);
m_1 = size(G_1,2);
m_2 = size(G_2,2);
q_1 = size(A_1,1);
q_2 = size(A_2,1);

% create convex hull ellipsotope
c_CH = 0.5*(c_1+c_2);
G_CH = [G_1, G_2, 0.5*(c_1-c_2), zeros(n,2*(m_1+m_2))];

A_31 = [eye(m_1); -eye(m_1); zeros(2*m_2,m_1)];

A_32 = [zeros(2*m_1,m_2) ; eye(m_2); -eye(m_2)];

A_30 = [-0.5*ones(2*m_1,1); 0.5*ones(2*m_2,1)];

A_CH = [A_1,            zeros(q_1,m_2), -b_1/2, zeros(q_1,2*(m_1+m_2));
        zeros(q_2,m_1), A_2,             b_2/2, zeros(q_1,2*(m_1+m_2));
        A_31,           A_32,            A_30,  eye(2*(m_1+m_2))];

b_CH = [b_1/2; b_2/2; -0.5*ones(2*(m_1+m_2),1)] ;

% J_extra = (m_1 + m_2 + 1):(3*(m_1 + m_2) + 1) ;
% I_extra = num2cell(J_extra) ;
% I_CH = [{1:m_1, m_1+1:m_1+m_2}, I_extra] ;
m_3 = m_1 + m_2 ;
J_extra = (m_3 + 2):(3*m_3 + 1) ;
I_extra = num2cell(J_extra) ;
I_CH = [{1:m_1}, {m_1+1:m_1+m_2}, {m_3+1}, {(1:m_1) + m_3 + 1}, {(m_1+1:m_1+m_2) + m_3 + 1}] ;

%% replicate with constrained zonotopes
CZ_1 = conZonotope(c_1,G_1,A_1,b_1) ;
CZ_2 = conZonotope(c_2,G_2,A_2,b_2) ;
CZ_CH = conZonotope(c_CH,G_CH,A_CH,b_CH) ;

%% plotting
figure(1); clf; axis equal; hold on; grid on;

% plot etopes
E1 = ellipsotope(p_norm,c_1,G_1,A_1,b_1);
E2 = ellipsotope(p_norm,c_2,G_2,A_2,b_2);
E_CH = ellipsotope(p_norm,c_CH,G_CH,A_CH,b_CH,I_CH);
plot(E1); plot(E2); plot(E_CH);

% figure(1) ; clf ; axis equal ; hold on ; grid on ;
% plot_zonotope({c_1,G_1,A_1,b_1})
% plot_zonotope({c_2,G_2,A_2,b_2})
% plot_zonotope({c_CH,G_CH,A_CH,b_CH},'facecolor','r','edgecolor','r')

%% helper functions
function h = plot_zonotope(zono_spec,varargin)
    % make zonotope with CORA
    [c,G,A,b] = get_zono_params_from_spec(zono_spec) ;
    Z = conZonotope(c,G,A,b) ;    
    
    % get vertices
    V = vertices(Z) ;
    
    if isempty(V)
        warning('Input zonotope is empty!')
    end
    
    % get faces
    F = 1:size(V,2) ;
    
    % plot with default input argz
    h = patch('faces',F,'vertices',V','facealpha',0.1,'facecolor','b','edgecolor','b',...
        varargin{:}) ;
    
    if nargout < 1
        clear h ;
    end
end

function [c,G,A,b] = get_zono_params_from_spec(zono_spec)
    c = zono_spec{1} ;
    G = zono_spec{2} ;
    if length(zono_spec) > 2
        A = zono_spec{3} ;
        b = zono_spec{4} ;
    else
        A = [] ;
        b = [] ;
    end
end