%% description
% This script demosntrates overapproximating an ellipsotope using
% [1, Prop. 5], which allows us to delete a constraint (wow!)
%
% [1] Scott, J.K., Raimondo, D.M., Marseglia, G.R. and Braatz, R.D., 2016.
%     Constrained zonotopes: A new tool for set-based estimation and fault
%     detection. Automatica, 69, pp.126-136.
%
% Authors: Shreyas Kousik
% Created: 29 May 2021
% Updated: 16 Mar 2022 (added flag to save figure, plus functionized stuff)
clear ; clc
%% user parameters
% rng
% rng(0)

% etope specs
p_norm = 2 ;
n_dim = 2 ;
n_gen = 8 ;
n_con = 2 ;

% whether or not to save output figure
flag_save_figure = false ;

%% autoamted from here
% make random etope
[E,c,G,A,b,I] = make_random_ellipsotope(p_norm,n_dim,n_gen,n_con) ;

% create Gamma and Lambda matrices
Gm = zeros(n_dim,n_con) ;
% Lm = 10*rand(n_con,n_con) - 5 ;
% Lm = diag([1,zeros(1,n_con-1)]) ;

%% iterate over the constraints and eliminate one at a time
E_list = cell(1,n_con) ;

for idx = 1:n_con
    E_list{idx} = reduce_etope_constraint(E,idx) ;
end

%% compare against reducing a constraint and a generator
E_rdc_con_and_gen = reduce_etope_constraint_and_generator(E) ;

%% plotting
fh = figure(1) ; clf ; axis equal ; hold on ; grid on ;

% plot all etopes
for idx = 1:n_con
    plot(E_list{idx},'color','r','linestyle','--','facealpha',0.1)
end

% plot con+gen reduced
plot(E_rdc_con_and_gen,'color','g','linestyle','--','facealpha',0.1)

% plot origetope
plot(E,'facecolor',[0.7 0.7 1],'facealpha',1) ;

% xlabel('x\langle1\rangle')
% ylabel('x\langle2\rangle')
set_plot_linewidths(2)
set_plot_fontsize(15)

%% save figure
if flag_save_figure
    save_figure_to_png(fh,'order_reduc_constraint_removal.png')
end