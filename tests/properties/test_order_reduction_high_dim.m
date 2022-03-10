%% description
% This script creates a figure illustrating the order reduction heuristic
% for 2-ellipsotopes in higher dimensions.
%
% For n_topes = 12 and n_dim = 100, this script takes about 13 s to run on
% a 2020 MacBook Pro.
%
% See also: figure_order_reduction_2_etopes.m
%
% Authors: Shreyas Kousik
% Created: 21 Feb 2022
% Updated: 5 Mar 2022
clear ; clc
%% user parameters
% random number generator seed
% rng(0) ;

% number of topes to test (things blow up quadratically so leave this less
% than, like, 20; default is 6)
n_comp_topes = 6 ; % default is 6

% dimension
n_dim = 10 ; % default is 10

%% automated from here
disp('Generating random ellipsotopes')

% create shared center
c = zeros(n_dim,1) ;

% % uncomment the following lines to make super-similar volume 'topes
% G = inv(sqrtm(make_random_covariance_matrix(n_dim))) ;
% G = G./max(G(:)) ;
% G = 2.*rand(n_dim) - 1 ;
 
% create all ellipsotopes
E_cell = cell(1,n_comp_topes) ;

for idx = 1:n_comp_topes
    % R_idx = RandOrthMat(n_dim) ; % make_random_orthonormal_matrix(n_dim) ;
    % s_idx = 1 ; %(2*rand(1) - 1) + 1 ;
    % G_idx = R_idx*G ;
    % s_idx = rand_range(0.1,1) ;
    
    % create generator matrix that is "well-scaled" (this prevents poor
    % conditioning of the ellipsoid shape matrix)
    G_idx = (1/sqrt(n_dim))*(2.*rand(n_dim) - 1) ;
    
    % create 'tope
    E_cell{idx} = ellipsotope(2,c,G_idx) ;
end

%% creating a new heuristic
disp('Computing true MVOE volumes and heuristic values')

% possible combinations of topes
combs = combinator(n_comp_topes,2,'c') ;
n_combs = size(combs,1) ;

% for each combination...
vols_MVOE = nan(1,n_combs) ; % volume of MVOE
vols_heur = nan(1,n_combs) ; % heuristic volume value
bts = nan(1,n_combs) ;

t_MVOE = nan(1,n_combs) ;
t_heur = nan(1,n_combs) ;

t_start = tic ;
for idx = 1:n_combs
    % get ellipsoids to compare
    idx_i = combs(idx,1) ;
    idx_j = combs(idx,2) ;
    
    % get the centers and generator matrices
    c_i = E_cell{idx_i}.center ;
    c_j = E_cell{idx_i}.center ;
    G_i = E_cell{idx_i}.generators ;
    G_j = E_cell{idx_j}.generators ;
    
    % compute the MVOE
    [G_rdc,~,bt] = make_MVOE_generator_matrix(G_i,G_j) ;
    vols_MVOE(idx) = ellipsoid_volume_from_generator_matrix(G_rdc) ;
    bts(idx) = bt ;
    
    % compute heuristic
    % get shape matrices AS IN THE HALDER PAPER! (inverse of etope paper)
    Q_i = inv(pinv(G_i)'*pinv(G_i)) ;
    Q_j = inv(pinv(G_j)'*pinv(G_j)) ;
    Q_MVOE = 2*(Q_i + Q_j) ; % assume bt = 1 ;
    vols_heur(idx) = 1/det(inv(Q_MVOE)) ;
    
    % measure timing
    t_heur(idx) = timeit(@() ellipsoid_volume_from_generator_matrix(G_rdc)) ;
    t_MVOE(idx) = timeit(@() make_MVOE_generator_matrix(G_i,G_j)) ;
end

toc(t_start)

%% test heuristic quality
disp('Testing heuristic quality')

% code from https://www.mathworks.com/help/matlab/data_analysis/linear-regression.html
x = vols_MVOE ;
y = vols_heur ;

p = polyfit(x,y,1) ;
yfit = polyval(p,x) ;
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal ;

%% plotting setup
disp('Plotting!')
x_ln = [0, max(vols_MVOE)] ;
y_ln = polyval(p,x_ln) ;

%% plotting
figure(1) ; clf ;

% MVOE vs. heuristic correlation
subplot(1,3,1) ; grid on ; axis equal ; hold on ;

plot(x_ln,y_ln,'-','Color',[1 0.5 0])
plot(vols_MVOE,vols_heur,'b.')

title('true vs. heuristic volume')
xlabel('true MVOE volume')
ylabel('heuristic value')
text(x_ln(1)+0.1,0.85*y_ln(2),['r^2 = ',num2str(rsq,'%0.4f')])
make_plot_pretty()

% zeta values
subplot(1,3,2)
histogram(bts)
title('MVOE \zeta values','interpreter','tex')
xlabel('bins')
ylabel('counts')
make_plot_pretty()

% timing
subplot(1,3,3) ; hold on ; grid on
boxplot([t_MVOE ; t_heur]','labels',{'MVOE','heuristic'},'labelorientation','horizontal')
n = findobj(gcf,'tag','Outliers');
for j = 1:numel(n)
    n(j).MarkerEdgeColor = 'b';
end
axis 'auto y'
ylabel('computation time [s]')
% set(gca, 'YScale', 'log')
make_plot_pretty()