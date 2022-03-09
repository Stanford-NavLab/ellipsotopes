%% description
% This script illustrates the order reduction heuristic for 2-ellipstopes
% over a variety of state dimensions and component ellipsoids.
%
% For 3 different dimensions, 6 component topes, and 100 topes per dim,
% this takes about 2 mins to run on a 2020 MacBook Pro.
%
% See also: figure_order_reduction_2_etopes.m,
% test_order_reduction_high_dim.m
%
% Authors: Shreyas Kousik
% Created: 5 Mar 2022
% Updated: --
clear ; clc ; close all ;
%% user parameters
% random number seed
rng(0)

% number of component ellipsoids
n_comp_topes = 6 ; % default is 10

% dimensions to test
n_dim_list = [2 8 14] ;

% number of topes to test
n_topes_per_dim = 100 ;

% whether or not to save figures
flag_save_figures = true ;

%% automated from here
disp('Testing MVOE volume vs. heuristic')

% possible combinations of topes
combs = combinator(n_comp_topes,2,'c') ;
n_combs = size(combs,1) ;
n_n_dim = length(n_dim_list) ;

% set up to save volume and timing info
vols_MVOE = nan(n_n_dim,n_combs) ; % volume of MVOE
vols_heur = nan(n_n_dim,n_combs) ; % heuristic volume value
bts = nan(n_n_dim,n_combs) ;

t_MVOE = nan(n_n_dim,n_combs) ;
t_heur = nan(n_n_dim,n_combs) ;

t_start = tic ;
for idx_dim = 1:n_n_dim
    % current dimension
    n_dim = n_dim_list(idx_dim) ;
    
    disp(['Testing n = ',num2str(n_dim)])
    
    for idx_tope = 1:n_topes_per_dim
        %% generate random etope from component etopes
        % create shared center
        c = zeros(n_dim,1) ;
        
        % create all ellipsotopes
        E_cell = cell(1,n_comp_topes) ;
        
        for idx_comp_topes = 1:n_comp_topes
            % create generator matrix that is "well-scaled" (this prevents
            % poot conditioning of the ellipsoid shape matrix)
            G_idx = (1/sqrt(n_dim))*(2.*rand(n_dim) - 1) ;
            
            % create component 'tope
            E_cell{idx_comp_topes} = ellipsotope(2,c,G_idx) ;
        end
        
        %% test MVOE for each combination of topes
        for idx_comb = 1:n_combs
            % get ellipsoids to compare
            idx_i = combs(idx_comb,1) ;
            idx_j = combs(idx_comb,2) ;
            
            % get the centers and generator matrices
            c_i = E_cell{idx_i}.center ;
            c_j = E_cell{idx_i}.center ;
            G_i = E_cell{idx_i}.generators ;
            G_j = E_cell{idx_j}.generators ;
            
            % compute the MVOE
            [G_rdc,~,bt] = make_MVOE_generator_matrix(G_i,G_j) ;
            vols_MVOE(idx_dim,idx_comb) = ellipsoid_volume_from_generator_matrix(G_rdc) ;
            bts(idx_dim,idx_comb) = bt ;
            
            % compute heuristic
            vols_heur(idx_dim,idx_comb) = compute_order_reduction_MVOE_heuristic(G_i,G_j) ;
            
            % measure timing
            t_heur(idx_dim,idx_comb) = timeit(@() compute_order_reduction_MVOE_heuristic(G_i,G_j)) ;
            t_MVOE(idx_dim,idx_comb) = timeit(@() make_MVOE_generator_matrix(G_i,G_j)) ;
        end
        
        % display progress
        idx_disp = round(100*idx_tope/n_topes_per_dim) ;
        if mod(idx_disp,10) == 0
            disp(['    ',num2str(idx_disp),' %'])
        end
    end
    
    disp(['    Time spent: ',num2str(toc(t_start),'%0.2f'),' s'])
end

% total time
t_total = toc(t_start) ;
disp(['Total time spent: ',num2str(t_total,'%0.2f'),' s'])

%% test heuristic quality
disp('Testing heuristic quality')

% set up to save
p_fit = nan(n_n_dim,2) ;
rsq_fit = nan(1,n_n_dim) ;

for idx = 1:n_n_dim
    % code from https://www.mathworks.com/help/matlab/data_analysis/linear-regression.html
    p = polyfit(vols_MVOE(idx,:),vols_heur(idx,:),1) ;
    y_fit = polyval(p,vols_MVOE(idx,:)) ;
    y_resid = vols_heur(idx,:) - y_fit;
    SS_resid = sum(y_resid.^2);
    SS_total = (length(vols_heur(idx,:))-1) * var(vols_heur(idx,:));
    rsq = 1 - SS_resid/SS_total ;
    
    % save data
    p_fit(idx,:) = p ;
    rsq_fit(idx) = rsq ;
end
%% plotting
disp('Plotting!')

for idx = 1:n_n_dim
    figure(idx) ; clf ;
    
    % create best fit line
    x_ln = [0, max(vols_MVOE(idx,:))] ;
    y_ln = polyval(p,x_ln) ;

    %% MVOE vs. heuristic correlation
    subplot(1,2,1) % subplot(1,3,1)
    grid on ; hold on ; axis square ;

    plot(x_ln,y_ln,'-','Color',[1 0.5 0])
    plot(vols_MVOE(idx,:),vols_heur(idx,:),'b.')

    xlabel('MVOE volume')
    ylabel('heuristic value')
    text(0.1*(x_ln(1)+x_ln(2)),0.85*y_ln(2),['r^2 = ',num2str(rsq,'%0.4f')])
    make_plot_pretty()

    %% timing
    subplot(1,2,2) ; % subplot(1,3,3)
    hold on ; grid on ; axis square ;
    boxplot([t_MVOE(idx,:) ; t_heur(idx,:)]','labels',{'MVOE','heuristic'},'labelorientation','horizontal')
    n = findobj(gcf,'tag','Outliers');
    for j = 1:numel(n)
        n(j).MarkerEdgeColor = 'b';
    end
    axis 'auto y'
    ylabel('computation time [s]')
    % set(gca, 'YScale', 'log')
    make_plot_pretty()
    
    
    %% zeta values
%     subplot(1,3,3) ; axis square ;
%     histogram(bts(idx,:))
%     xlabel('MVOE \zeta values','interpreter','tex')
%     ylabel('counts')
%     make_plot_pretty()

    %% save figure
    if flag_save_figures
        filename = ['order_reduction_MVOE_heuristic_',num2str(n_dim_list(idx)),'D.png'] ;
        exportgraphics(gcf,filename)
    end
end