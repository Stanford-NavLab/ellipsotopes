%% user parameters
% p norm
p = 6 ;

%% automated from here
% make a circle of points
[F,V] = make_circle(1,100) ;

% given the p-norm, make a p-ircle
x_1_vec = linspace(-1,1,50) ;
x_2_vec_up = (1 - x_1_vec.^p).^(1/p) ;
x_2_vec_dn = -x_2_vec_up ;
E = [x_1_vec, x_1_vec(end-1:-1:1) ;
     x_2_vec_up, x_2_vec_dn(end-1:-1:1)] ;

%% make "ground truth" for comparison
% make points
X = make_grid_2D() ;

% evaluate p-norm on points
p_norm_log = sum(X.^p,1) <= 1 ;
X = X(:,p_norm_log) ;
K = boundary(X') ;
X = X(:,K) ;

%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ;

plot_path(X,'b--')
plot_path(E,'r-')

set_plot_linewidths(1.5)