%% description
% Test plotting several example ellipsotopes using the built-in class
% method.
%
% Authors: Adam Dai (Shreyas is sneakin in here too)
% Created: 14 Apr 2021
% Updated: 15 Apr 2021 (added fontsize settings to each plot lol)
clear ; clc ;
%% user parameters
% rng seed (for random etopes)
rng(0)

% set to true to save figures
flag_save_figures_to_png = true ;

%% basic 2D 2-ellipsotope
p = 2;
c = zeros(2,1);
G = [2 1; 1 2];

fh = figure(); axis equal; grid on
E = ellipsotope(p,c,G);
plot(E); title('Basic 2D 2-ellipsotope');
set(gca,'fontsize',15)
if flag_save_figures_to_png
    save_figure_to_png(fh,'plot_ex_basic_2D_2_ellipsotope.png')
end

%% constrained 2D 2-ellipsotope
p = 2;
c = zeros(2,1);
G = [1 -1 0.2; 1 2 1];
A = [1 -1 1];
b = 0.5;

fh = figure(); axis equal; grid on
E = ellipsotope(p,c,G,A,b);
plot(E); title('Constrained 2D 2-ellipsotope');
set(gca,'fontsize',15)
if flag_save_figures_to_png
    save_figure_to_png(fh,'plot_ex_constrained_2D_2_ellipsotope.png')
end

%% indexed 2D 2-ellipsotope
p = 2;
c = zeros(2,1);
G = [1 -1 0.2; 1 2 1];
I = {[1 2],3};

fh = figure(); axis equal; grid on
E = ellipsotope(p,c,G,[],[],I);
plot(E); title('Indexed 2D 2-ellipsotope');
set(gca,'fontsize',15)
if flag_save_figures_to_png
    save_figure_to_png(fh,'plot_ex_indexed_2D_2_ellipsotope.png')
end

%% zonotope using indexed ellipsotope (eqn 8)
p = 2;
c = zeros(2,1);
G = [1 -1 0.2; 1 2 1];
I = {1,2,3};

fh = figure(); axis equal; grid on
E = ellipsotope(p,c,G,[],[],I);
plot(E); title('Zonotope using indexed ellipsotope (eqn 8)')
set(gca,'fontsize',15)
if flag_save_figures_to_png
    save_figure_to_png(fh,'plot_ex_zonotope_as_ellipsotope.png')
end

%% random 2D general 2-ellipsotope
p = 2;
c = zeros(2,1);
G = 2*rand(2,4) - 1;
A = [-1 1 -1 1];
b = 0.5;
I = {[1,2],[3,4]};

fh = figure(); axis equal; grid on
E = ellipsotope(p,c,G,A,b,I);
plot(E); title('Random 2D general ellipsotope');
set(gca,'fontsize',15)
if flag_save_figures_to_png
    save_figure_to_png(fh,'plot_ex_random_ellipsotope.png')
end

%% basic 2D 3-ellipsotope
p = 3;
c = zeros(2,1);
G = [2 1; 1 2];

fh = figure(); axis equal; grid on
E = ellipsotope(p,c,G);
plot(E); title('Basic 2D 3-ellipsotope');
set(gca,'fontsize',15)

%% basic 3D 2-ellipsotope
p = 2;
c = zeros(3,1);
G = eye(3);

fh = figure(); axis equal; grid on; 
E = ellipsotope(p,c,G);
plot(E); title('Basic 3D 2-ellipsotope');
set(gca,'fontsize',15)
if flag_save_figures_to_png
    save_figure_to_png(fh,'plot_ex_basic_3D_2_ellipsotope.png')
end

%% general 3D 2-ellipsotope
p = 2;
c = zeros(3,1);
G = [1 0 0 0.5;
     0 1 0 0;
     0 0 1 0.5];
A = [-1 1 -1 0.5];
b = 0.5;
I = {[1,2],[3 4]};

fh = figure(); axis equal; grid on; 
E = ellipsotope(p,c,G,A,b,I);
plot(E,'facealpha',1); shading interp; title('General 3D 2-ellipsotope');
set(gca,'fontsize',15)
if flag_save_figures_to_png
    save_figure_to_png(fh,'plot_ex_3D_2_ellipsotope.png')
end

%% x-y projection of 3D 2-ellipsotope
p = 2;
c = zeros(3,1);
G = [1 0 0 0.5;
     0 1 0 0;
     0 0 1 0.5];
A = [-1 1 -1 0.5];
b = 0.5;
I = {[1,2],[3 4]};

fh = figure(); axis equal; grid on; 
E = ellipsotope(p,c,G,A,b,I);
%E = ellipsotope(p,c,G);
plot(E,'proj_dims',[1 2]); title('x-y projection of 2-ellipsotope');
set(gca,'fontsize',15)
if flag_save_figures_to_png
    save_figure_to_png(fh,'plot_ex_3D_2_ellipsotope_projection.png')
end