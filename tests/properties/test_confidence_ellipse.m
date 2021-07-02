%% description
% This script tests generates the confidence ellipse (defined by a 
% probability threshold) of a Gaussian.
%
% Authors: Adam Dai 
% Created: 20 May 2021
% Updated: 
%
%% user parameters
% Gaussian
mu = [0;0];
%Sigma = [2 1;1 2];
% random Sigma
A = rand(2,2);
Sigma = A'*A;

% for 2D
P = 0.9973; % probability threshold
eps = -2*log(1-P);
G = sqrtm(eps*Sigma);
E = ellipsotope(2,mu,G);

% whether or not to save flag
flag_save_figure = false;

%% automated from here

% sample points from the distribution
n_samples = 1e4;
P = mvnrnd(mu,Sigma,n_samples);
 
% test how many points fall in the confidence ellipse
n_inside = 0;
for i = 1:n_samples
    p = P(i,:)';
    if E.contains(p)
        n_inside = n_inside + 1;
    end
end

disp(['N(inside):', num2str(n_inside)]);
disp(['P(inside):', num2str(n_inside/n_samples)]);

%% plotting
h = figure(1) ; clf ; axis equal ; hold on ; grid on ;
plot(E,'facecolor','r','edgecolor','r','facealpha',0.1);
scatter(P(:,1),P(:,2),'.');
lim = axis; axis(lim + [-1 1 -1 1]);

if flag_save_figure
    save_figure_to_pdf(h,'confidence_ellipse.pdf')
    save_figure_to_png(h,'confidence_ellipse.png')
end