%% description
% This script evaluates the emptiness check for ellipsotopes as a function
% of the number of generators. We check the proposition versus the
% corollary to see which is faster.
%
% NOTE this takes about 25 s to run on an i7 8-core processor.
%
% Authors: Shreyas Kousik
% Created: 27 Jul 2021
% Updated: -- 
clear ; clc
%% user parameters
% rng seed
rng(0)

% properties that stay fixed
p_norm = 2 ;
n_dim = 2 ;
n_con = 1 ;

% range of properties to test
n_gen_list = 1:10:101 ; % default is 1:20
n_etopes_per_n_gen = 5 ; % default is 10

%% automated from here
n_n_gen = length(n_gen_list) ;

% set up to store times
t_avg_prop = nan(n_etopes_per_n_gen,n_n_gen) ; % for nonempty etopes
t_avg_empty = nan(n_etopes_per_n_gen,n_n_gen) ; % for empty etopes

%% time etope emptiness check
disp('Timing emptiness checks!')
start_tic_all = tic ;
for idx_gen = 1:n_n_gen 
    n_gen = n_gen_list(idx_gen) ;
    disp(['Testing ',num2str(n_gen),' generators']) ;
    start_tic_gen = tic ;
    for idx_tope = 1:n_etopes_per_n_gen
        % make a random ellipsotope
        [E,c,G,A,~,I] = make_random_ellipsotope(p_norm,n_dim,n_gen,n_con) ;
        E.constraint_b = 10*n_gen.*rand(n_con,1) ;

        %% run emptiness checks
        t_avg_prop(idx_tope,idx_gen) = timeit(@() E.isempty(true,'standard')) ;
        t_avg_empty(idx_tope,idx_gen) = timeit(@() E.isempty(true,'feasibility')) ;
    end
    toc(start_tic_gen)
end

disp([newline,'Total time elapsed:',newline,num2str(toc(start_tic_all),'%0.1f'),' s']) ;

%% plotting
fh = figure(1) ; clf ;

subplot(2,1,1) ; hold on ; grid on ;
h_prop = boxplot(t_avg_prop,'Colors','b','labels',n_gen_list,'labelorientation','horizontal') ;
axis 'auto y'

n = findobj(gcf,'tag','Outliers');
for j = 1:numel(n)
    n(j).MarkerEdgeColor = 'b';
end

xlabel('# of generators')
ylabel('time [s]')
set_plot_linewidths(1.5) ;
set_plot_fontsize(15) ;

subplot(2,1,2) ; hold on ; grid on ;
h_cor = boxplot(t_avg_empty,'Colors','r','labels',n_gen_list,'labelorientation','horizontal') ;
axis 'auto y'

% legend([h_full,h_empty],{'full','empty'},'location','northwest')
xlabel('# of generators')
ylabel('time [s]')
set_plot_linewidths(1.5) ;
set_plot_fontsize(15) ;

save_figure_to_png(fh,'emptiness_check_time_prop_vs_cor.png')