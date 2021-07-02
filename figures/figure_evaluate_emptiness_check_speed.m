%% description
% This script evaluates the emptiness check for ellipsotopes as a function
% of the number of generators.
%
% NOTE this takes about 25 s to run on an i7 8-core processor.
%
% Authors: Shreyas Kousik
% Created: 28 May 2021
% Updated: 15 June 2021 (reran with bug fixed in E.isempty())
clear ; clc
%% user parameters
% rng seed
rng(0)

% properties that stay fixed
p_norm = 2 ;
n_dim = 2 ;
n_con = 1 ;

% range of properties to test
n_gen_list = 1:20 ; % default is 1:20
n_etopes_per_n_gen = 10 ; % default is 10

%% automated from here
n_n_gen = length(n_gen_list) ;

% set up to store times
t_avg_full = nan(n_etopes_per_n_gen,n_n_gen) ; % for nonempty etopes
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

        %% time for nonempty etope
        % set b for nonempty etopes
        E.constraint_b = zeros(n_con,1) ;

        % run emptiness check
        t_avg = timeit(@() E.isempty()) ;
        
        t_avg_full(idx_tope,idx_gen) = t_avg ;
        
        %% time for empty etope
        % set b for empty etopes
        E.constraint_b = 10.*ones(n_con,1) ;

        % run emptiness check
        t_avg = timeit(@() E.isempty()) ;
        
        t_avg_empty(idx_tope,idx_gen) = t_avg ;
    end
    toc(start_tic_gen)
end

disp([newline,'Total time elapsed:',newline,num2str(toc(start_tic_all),'%0.1f'),' s']) ;

%% plotting
fh = figure(1) ; clf ;

subplot(2,1,1) ; hold on ; grid on ;
h_full = boxplot(t_avg_full,'Colors','b','labels',n_gen_list,'labelorientation','horizontal') ;
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
h_empty = boxplot(t_avg_empty,'Colors','r','labels',n_gen_list,'labelorientation','horizontal') ;
axis 'auto y'

% legend([h_full,h_empty],{'full','empty'},'location','northwest')
xlabel('# of generators')
ylabel('time [s]')
set_plot_linewidths(1.5) ;
set_plot_fontsize(15) ;

save_figure_to_png(fh,'emptiness_check_time.png')