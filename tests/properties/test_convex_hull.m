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
E1 = ellipsotope(p_norm,c_1,G_1,A_1,b_1);
E2 = ellipsotope(p_norm,c_2,G_2,A_2,b_2);

E_CH = conv_hull(E1,E2);

%% plotting
figure(1); clf; axis equal; hold on; grid on;

% plot etopes
plot(E1); plot(E2); plot(E_CH);

% %% replicate with constrained zonotopes
% CZ_1 = conZonotope(c_1,G_1,A_1,b_1) ;
% CZ_2 = conZonotope(c_2,G_2,A_2,b_2) ;
% CZ_CH = conZonotope(c_CH,G_CH,A_CH,b_CH) ;
% 
% % figure(1) ; clf ; axis equal ; hold on ; grid on ;
% % plot_zonotope({c_1,G_1,A_1,b_1})
% % plot_zonotope({c_2,G_2,A_2,b_2})
% % plot_zonotope({c_CH,G_CH,A_CH,b_CH},'facecolor','r','edgecolor','r')
% 
% %% helper functions
% function h = plot_zonotope(zono_spec,varargin)
%     % make zonotope with CORA
%     [c,G,A,b] = get_zono_params_from_spec(zono_spec) ;
%     Z = conZonotope(c,G,A,b) ;    
%     
%     % get vertices
%     V = vertices(Z) ;
%     
%     if isempty(V)
%         warning('Input zonotope is empty!')
%     end
%     
%     % get faces
%     F = 1:size(V,2) ;
%     
%     % plot with default input argz
%     h = patch('faces',F,'vertices',V','facealpha',0.1,'facecolor','b','edgecolor','b',...
%         varargin{:}) ;
%     
%     if nargout < 1
%         clear h ;
%     end
% end
% 
% function [c,G,A,b] = get_zono_params_from_spec(zono_spec)
%     c = zono_spec{1} ;
%     G = zono_spec{2} ;
%     if length(zono_spec) > 2
%         A = zono_spec{3} ;
%         b = zono_spec{4} ;
%     else
%         A = [] ;
%         b = [] ;
%     end
% end