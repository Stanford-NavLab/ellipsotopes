%% description
% This script tests solving for a scalar to push points to the boundary of
% the intersection of a linear subspace with a ball product.
%
% Authors: Shreyas Kousik
% Created: 6 Apr 2021
% Updated: 9 Apr 2021 (functionized stuff!)

clear;clc
%% user parameters
% rng seed
rng(0) ;

% p norm
p_norm = 2 ;

% number of points to create
n_P = 500 ;

% index set
I = {[1,2],[3]} ;

% constraint (it's prettier with just one constraint, but the code below
% works for an arbitrary number of constraints)
% A = [1 1 1 ;
%     -1 -1 1] ;
% b = 0.1*ones(2,1) ;
A = [1 -1 1] ;
b = 0.1 ;

% plotting projection dimensions
proj_dims = 1:3 ; 

%% automated from here
% TO DO: functionize these sanity checks within the ellipsotope.m class as
% a method

% sanity check the index set
[I_chk,n_I,n_dim] = check_index_set_validity(I) ;
if ~I_chk
    error(['The index set is not valid! ',...
        'It should contain each dimension of the hyperball product ',...
        'exactly once.'])
end

% get number of constraints
[n_con,n_dim_con] = size(A) ;

% sanity check the constraints
if n_dim_con ~= n_dim
    error(['The constraints defined by (A,b) must have as many columns ',...
        'as there are dimensions in the index set (i.e., the hyperball ',...
        'space of the ellipsotope.'])
end

%% project points to intersection of constraint and ball product
% create some random points
P_in = 2*rand(n_dim,n_P) - 1 ;

start_tic = tic ;
P_out = project_points_to_ball_product_and_linear_subspace(P_in,p_norm,A,b,I) ;
total_time = toc(start_tic) ;

disp(['Time spent projecting points: ',num2str(total_time),' s'])

%% plotting setup
% create ball product for patch
[F_E,V_E] = make_ball_product_for_patch(p_norm,I,proj_dims) ;

% get nullspace and translation of each linear constraint
K_each = cell(1,n_con) ;
t_each = cell(1,n_con) ;
for idx = 1:n_con
    K_each{idx} = null(A(idx,:)) ;
    t_each{idx} = A(idx,:)\b(idx) ;
end

%% plotting
fh = figure(1) ; clf ; axis equal ; grid on ; hold on ; view(3) ;

% plot ball product
h_E = patch('faces',F_E,'vertices',V_E,'facealpha',0.1','edgealpha',0,'facecolor','b') ;

% % plot t_con
% h_t = plot_path(t_con,'b^','markerfacecolor','b') ;

% plot points on boundary
h_P = plot_path(P_out(proj_dims,:),'b.','markersize',8) ;

% plot hyperplanes
if n_dim == 3
    h_H = plot_planes_3D(A,b,2.5) ;
    legend([h_E,h_H(1),h_P],{'ball product','linear subspace','boundary of intersection'})
else
    legend([h_E, h_P],{'ball product','boundary of intersection'})
end

% plotting cleanup
xlabel('x_1')
ylabel('x_2')
zlabel('x_3')
set(gca,'fontsize',15) ;

% save_figure_to_png(fh,'ellipsotope_ball_product_int_lin_space_bdry.png')