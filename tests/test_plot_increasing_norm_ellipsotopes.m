%% description
% This script tests plotting ellipsotopes with more than 3 generators
%
% Authors: Shreyas Kousik
% Created: 15 Mar 2021
% Updated: 17 Mar 2021
%
%% user parameters
% list of norms
p_list = 2:2:10 ;

% rng seed
rng(10) ;

% center
c = rand(2,1) ;

% generators
G = 2*rand(2,10) - 1 ;

%% automated from here
% list of ellipsotopes
E = cell(1,length(p_list)) ;

% make ellipsotopes
idx = 1 ;
for p = p_list
    E{idx} = ellipsotope(p,c,G) ;
    idx = idx + 1 ;
end

% make comparable zonotope (requires CORA)
Z = zonotope(c,G) ;

%% plotting
h = figure(1) ; clf ; axis equal ; hold on ; grid on ;

plot(Z) ;
for idx = length(p_list):-1:1
    plot(E{idx}) ;
    pause(0.05)
end

xlabel('x_1')
ylabel('x_2')

set(gca,'fontsize',15)
% save_figure_to_png(h,'ellipsotope_increasing_norm.png') ;