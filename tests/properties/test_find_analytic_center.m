%% description
% This script computes the analytic center of an ellipsotope using the
% method in:
%
% S. Boyd and L. Vandenberghe, Convex optimization. Cambridge university press, 2004.
%
% Authors: Shreyas Kousik
% Created: 23 Dec 2021
% Updated: nah
clear ; clc ;
%% user parameters
% % rng seed
% rng(1)

% ellipsotope params
p_norm = 2 ;
n_dim = 2 ;
n_gen = 12 ;
n_con = 4 ;

%% automated from here
% make etope
[E,c,G,A,b,I] = make_random_ellipsotope(p_norm,n_dim,n_gen,n_con) ;

%% compute analytic center
% set up optimization options
options = optimoptions('fmincon','Display','off',...
    'SpecifyObjectiveGradient',true) ;

% solve for center with fmincon
x_0 = pinv(A)*b ; % initial guess
x = fmincon(@(x) cost_fn(x,p_norm,I),x_0,[],[],A,b,[],[],[],options) ;

% get center in workspace
p = c + G*x ;

%% plotting
if n_dim == 2
    figure(1) ; clf ; axis equal ; hold on ; grid on ;
    
    plot(E)
    plot(p(1),p(2),'rp')
    
    make_plot_pretty()
end

%% helper functions
function [c,gc] = cost_fn(x,p_norm,I)
    n_X = length(x) ;
    n_I = length(I) ;
    h = nan(n_I,1) ;
    gc = nan(1,n_X) ;
    for idx = 1:n_I
        J = I{idx} ;
        h(idx) = sum(x(J).^p_norm) - 1 ;
        gc(J) = -(p_norm*x(J).^(p_norm-1))./h(idx) ;
    end

    c = -sum(log(-h)) ;
end