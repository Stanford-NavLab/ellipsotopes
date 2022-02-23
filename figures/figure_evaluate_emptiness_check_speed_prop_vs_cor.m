%% description
% This script evaluates the emptiness check for ellipsotopes as a function
% of the number of generators. We check the proposition versus the
% corollary to see which is faster.
%
% NOTE this takes ~52 s to run on an i7 8-core processor with the default
% user parameters (see comments below). The runtime is a lot more for
% any p > 2, which makes sense since we have to switch from solving a QP to
% solving a generic nonlinear (convex) program.
%
% We strongly recommend not bumping up the number of generators much, as
% the emptiness check can take a long time to solve for an empty
% ellipsotope with more than 20 generators.
%
% Authors: Shreyas Kousik
% Created: 27 Jul 2021
% Updated: 23 Feb 2022 (fiddled with plotting and stuff)
clear ; clc
%% user parameters
% rng seed
rng(0)

% properties that stay fixed
p_norm = 2 ; % default is 2
n_con = 1 ; % default is 1

% range of properties to test
n_dim_list = [2 8 14] ; % default is 2:2:10
n_gen_list = 2:2:20 ; % default is 2:2:10
n_etopes_per_n_gen = 5 ; % default is 5

% whether or not to save the final figure
flag_save_figure = true ;

%% automated from here
n_n_dim = length(n_dim_list) ;
n_n_gen = length(n_gen_list) ;

% set up to store times
t_avg_prop = nan(n_etopes_per_n_gen,n_n_gen,n_n_dim) ; % for nonempty etopes
t_avg_cor = nan(n_etopes_per_n_gen,n_n_gen,n_n_dim) ; % for empty etopes

%% time etope emptiness check
disp(['Timing emptiness checks!',newline])

start_tic_all = tic ;
parfor idx_dim = 1:n_n_dim
    n_dim = n_dim_list(idx_dim) ;
    disp([newline,'Testing dimension ',num2str(n_dim)]) ;
    
    for idx_gen = 1:n_n_gen 
        n_gen = n_gen_list(idx_gen) ;
        disp(['   Testing ',num2str(n_gen),' generators']) ;
        start_tic_gen = tic ;
        
        for idx_tope = 1:n_etopes_per_n_gen
            % make a random ellipsotope
            [E,c,G,A,~,I] = make_random_ellipsotope(p_norm,n_dim,n_gen,n_con) ;
            E.constraint_b = 10*n_gen.*rand(n_con,1) ;

            %% run emptiness checks
            t_avg_prop(idx_tope,idx_gen,idx_dim) = timeit(@() E.isempty(true,'standard')) ;
            t_avg_cor(idx_tope,idx_gen,idx_dim) = timeit(@() E.isempty(true,'feasibility')) ;
        end
        t_spent = toc(start_tic_gen) ;
        disp(['       Time elapsed: ',num2str(t_spent,'%0.3f'),' s'])
    end
end

disp([newline,'Total time elapsed:',newline,num2str(toc(start_tic_all),'%0.1f'),' s']) ;

%% plotting
% set up colors for the different dimensions
t_vec = linspace(0,pi,n_n_dim) ;
colors = 0.4.*[sin(t_vec); sin(t_vec + pi/5); sin(t_vec + 4*pi/3)]' + 0.4 ;
d_dim_offset = 0.2 ;
d_offset = (n_per_n_gen/5) ;
n_per_n_gen = mean(diff(n_gen_list)) ;

% start plottin'
fh = figure(1) ; clf ;

subplot(2,1,1) ; hold on ; grid on ;

for idx = 1:n_n_dim
    m = mean(t_avg_prop(:,:,idx),1) ;
    m_min = m - min(t_avg_prop(:,:,idx),[],1) ;
    m_max = m + max(t_avg_prop(:,:,idx),[],1) ;
    errorbar(n_gen_list - d_offset + d_dim_offset.*(idx),m,m_min,m_max,'.',...
        'Color',colors(idx,:),'linewidth',2,'CapSize',0,'MarkerSize',20)
end
legend(arrayfun(@(n) ['n = ',num2str(n)],n_dim_list,'Uni',0),'location','northwest')
xlabel('# of generators')
ylabel('Prop. 7 time [s]')

make_plot_pretty()

subplot(2,1,2) ; hold on ; grid on ;
for idx = 1:n_n_dim
    m = mean(t_avg_cor(:,:,idx),1) ;
    m_min = m - min(t_avg_cor(:,:,idx),[],1) ;
    m_max = m + max(t_avg_cor(:,:,idx),[],1) ;
    errorbar(n_gen_list - d_offset + d_dim_offset.*(idx),m,m_min,m_max,'.',...
        'Color',colors(idx,:),'linewidth',2,'CapSize',0,'MarkerSize',20)
end
xlabel('# of generators')
ylabel('Cor. 8 time [s]')

make_plot_pretty()

if flag_save_figure
    save_figure_to_png(fh,'emptiness_check_time_prop_vs_cor.png')
end