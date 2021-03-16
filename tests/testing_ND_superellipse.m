%% description
% This script plots points sampled randomly from an n-dimensional
% superellipse. Since we can't plot n dimensions, we plot a bunch of 2-D
% projections of the points.
%
% Authors: Shreyas Kousik
% Created: 15 Mar 2021
% Updated: not yet!
%
%% user parameters
p = 8 ; % norm
d = 5 ; % dimension
n = 10000 ; % number of points

%% automated from here
% generate n points
P = 2.*rand(d,n) - 1 ;

% map all points to the p-sphere
C = P./repmat(vecnorm(P,p),d,1) ;

%% plotting
figure(1) ; clf ;

n_plot = ceil(sqrt(d-1)) ;

for idx = 1:(d-1)
    subplot(n_plot,n_plot,idx) ;
    axis equal ; hold on ;
    
    plot_path(P(idx:idx+1,:),'r.')
    plot_path(C(idx:idx+1,:),'b.')
    
    xlabel(['x_',num2str(idx)])
    ylabel(['x_',num2str(idx+1)])
    
    set(gca,'fontsize',12)
end