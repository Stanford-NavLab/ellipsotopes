%% description
% This script tests plotting ellipsotopes with more than 3 generators
%
% Authors: Shreyas Kousik
% Created: 15 Mar 2021
% Updated: 9 July 2021 (got ride of axis labels)
%
%% user parameters
% list of norms
p_list = [2:2:10, inf] ;

% rng seed
rng(1) ;

% center
c = rand(2,1) ;

% generators
n_gen = 5 ;
G = 2*rand(2,n_gen) - 1 ;

% rotate G to be horizontal
G = rotation_matrix_2D(-pi/6) * G; 

%% automated from here
% list of ellipsotopes
E = cell(1,length(p_list)) ;

% make ellipsotopes
idx = 1 ;
for p = p_list
    E{idx} = ellipsotope(p,c,G) ;
    idx = idx + 1 ;
end

%% plotting
h = figure(1) ; clf ; hold on ; grid on ; axis equal ; axis tight

for idx = length(p_list):-1:1
    plot(E{idx},'plot_method','ray','edgealpha',1.0,'linewidth',1) ;
    pause(0.05)
end

% xlabel('$x\langle 1 \rangle$','Interpreter','latex'); ylabel('$x\langle 2 \rangle$','Interpreter','latex');
%set(gca,'fontsize',15)
set_plot_fontsize(15)
lim = axis; axis(lim + 0.25*[-1 1 -1 1]);

%% saving
save_figure_to_png(h,'ellipsotope_increasing_norm.png') ;
save_figure_to_pdf(h,'ellipsotope_increasing_norm.pdf') ;