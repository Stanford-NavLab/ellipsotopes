%% description
% This script evaluates the emptiness check for ellipsotopes as a function
% of the number of generators.
%
% NOTE this takes about 25 s to run on an i7 8-core processor.
%
% Authors: Shreyas Kousik
% Created: 28 May 2021
% Updated: 23 Feb 2022 (updated plotting)
clear ; clc
%% user parameters
% rng seed
rng(0)

% properties that stay fixed
p_norm = 2 ;
n_con = 1 ;

% range of properties to test
n_dim_list = [2 8 14]; % default is 2:2:16
n_gen_list = 1:20:101 ; % 1:20:101 ; % default is 1:20:101
n_etopes_per_n_gen = 8 ; % default is 10

% whether or not to save the final figure
flag_save_figure = false ;

%% automated from here
n_n_dim = length(n_dim_list) ;
n_n_gen = length(n_gen_list) ;

% set up to store times
t_avg_full = nan(n_etopes_per_n_gen,n_n_gen,n_n_dim) ; % for nonempty etopes
t_avg_empty = nan(n_etopes_per_n_gen,n_n_gen,n_n_dim) ; % for empty etopes

%% time etope emptiness check
disp('Timing emptiness checks!')
start_tic_all = tic ;
for idx_dim = 1:n_n_dim
    n_dim = n_dim_list(idx_dim) ;
    disp([newline,'Testing dimension ',num2str(n_dim)]) ;
    
    for idx_gen = 1:n_n_gen
        n_gen = n_gen_list(idx_gen) ;
        disp([newline,'Testing ',num2str(n_gen),' generators']) ;
        start_tic_gen = tic ;
        for idx_tope = 1:n_etopes_per_n_gen
            % make a random ellipsotope
            [E,c,G,A,~,I] = make_random_ellipsotope(p_norm,n_dim,n_gen,n_con) ;
            
            %% time for nonempty etope
            % set b for nonempty etopes
            E.constraint_b = zeros(n_con,1) ;
            
            % run emptiness check
            t_avg = timeit(@() E.isempty(false,'feasibility')) ;
            
            t_avg_full(idx_tope,idx_gen,idx_dim) = t_avg ;
            
            %% time for empty etope
            % set b for empty etopes
            E.constraint_b = (2*n_gen).*ones(n_con,1) ;
            
            % % check that it's empty, yikes
            % assert(E.isempty(),'Something went wrong! Ellipsotope is nonempty')
            
            % run emptiness check
            t_avg = timeit(@() E.isempty(false,'feasibility')) ;
            
            t_avg_empty(idx_tope,idx_gen,idx_dim) = t_avg ;
        end
        toc(start_tic_gen)
    end
end

disp([newline,'Total time elapsed:',newline,num2str(toc(start_tic_all),'%0.1f'),' s']) ;

%% plotting
% set up colors for the different dimensions
t_vec = linspace(0,pi,n_n_dim) ;
colors = 0.4.*[sin(t_vec); sin(t_vec + pi/5); sin(t_vec + 4*pi/3)]' + 0.4 ;
d_dim_offset = 2 ;
d_offset = (n_per_n_gen/2) ;
n_per_n_gen = mean(diff(n_gen_list)) ;

fh = figure(1) ; clf ;

% subplot(2,1,1) ; hold on ; grid on ;
% h_full = boxplot(t_avg_full,'Colors','b','labels',n_gen_list,'labelorientation','horizontal') ;
% axis 'auto y'
% 
% n = findobj(gcf,'tag','Outliers');
% for j = 1:numel(n)
%     n(j).MarkerEdgeColor = 'b';
% end


subplot(2,1,1) ; hold on ; grid on ;

for idx = 1:n_n_dim
    % h_prop = boxplot(t_avg_prop(:,:,idx),'Colors',rand(1,3),'labels',n_gen_list,'labelorientation','horizontal') ;
    m = mean(t_avg_full(:,:,idx),1) ;
    m_min = m - min(t_avg_full(:,:,idx),[],1) ;
    m_max = m + max(t_avg_full(:,:,idx),[],1) ;
    errorbar(n_gen_list - (n_per_n_gen/4) + d_dim_offset.*(idx),m,m_min,m_max,'.',...
        'Color',colors(idx,:),'linewidth',2,'CapSize',0,'MarkerSize',20)
end
legend(arrayfun(@(n) ['n = ',num2str(n)],n_dim_list,'Uni',0),'location','northwest')

xlabel('# of generators')
ylabel('time [s]')

make_plot_pretty()
% 
% subplot(2,1,2) ; hold on ; grid on ;
% h_empty = boxplot(t_avg_empty,'Colors','r','labels',n_gen_list,'labelorientation','horizontal') ;
% axis 'auto y'
% 
% % legend([h_full,h_empty],{'full','empty'},'location','northwest')
% xlabel('# of generators')
% ylabel('time [s]')
% set_plot_linewidths(1.5) ;
% set_plot_fontsize(15) ;

subplot(2,1,2) ; hold on ; grid on ;

for idx = 1:n_n_dim
    % h_prop = boxplot(t_avg_prop(:,:,idx),'Colors',rand(1,3),'labels',n_gen_list,'labelorientation','horizontal') ;
    m = mean(t_avg_empty(:,:,idx),1) ;
    m_min = m - min(t_avg_empty(:,:,idx),[],1) ;
    m_max = m + max(t_avg_empty(:,:,idx),[],1) ;
    errorbar(n_gen_list - (n_per_n_gen/4) + d_dim_offset.*(idx),m,m_min,m_max,'.',...
        'Color',colors(idx,:),'linewidth',2,'CapSize',0,'MarkerSize',20)
end
legend(arrayfun(@(n) ['n = ',num2str(n)],n_dim_list,'Uni',0),'location','northwest')

xlabel('# of generators')
ylabel('time [s]')

make_plot_pretty()


save_figure_to_png(fh,'emptiness_check_time.png')