%% description
% This script tests the ellipsotope point containment lemma.
%
% Authors: Shreyas Kousik
% Created: 15 Apr 2021 (just barely, it's like 11:49 pm lol)
% Updated: not yet
clear ; clc
%% user parameters
% rng seed
rng(0)

% point to test
p_test = [1.0 ; 1.0] ; % this point is inside
% p_test = [1.5; 0.0] ; % this point is inside


% ellipsotope
p_norm = 4 ;
c = zeros(2,1) ;
G = rand(2,3) ;
A = [] ;
b = [] ;
I = {1:2,3} ;

%% automated from here
% get sizes of things
n_gen = size(G,2) ;
n_con = size(A,1) ;

% compute point containtment constraint matrices
A_eq = [G ; A] ;
b_eq = [p_test - c ; b] ;

% set up program
cost = @(x) point_cointainment_cost(x,p_norm,I) ;
x_0 = rand(n_gen,1) ;
options = optimoptions('fmincon','Display','off','CheckGradients',true) ;
[x_opt,f_val] = fmincon(cost,x_0,[],[],A_eq,b_eq,[],[],[],options) ;

if f_val <= 1
    disp('Test point is INSIDE the ellipsotope!')
else
    disp('Test point is OUTSIDE the ellipsotope!')
end

% time the containment check
avg_time = timeit(@() fmincon(cost,x_0,[],[],A_eq,b_eq,[],[],[],options)) ;
disp(['It takes ',num2str(avg_time,'%0.4f'),' s to compute the point containment.'])

%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ;

% plot etope
E = ellipsotope(p_norm,c,G,A,b,I) ;
plot(E) ;

% plot point
plot_path(p_test,'rx','linewidth',1.5)

% cleanup
legend('ellipsotope','test point','location','northwest')
set(gca,'fontsize',15)

%% helper functions
function [c,gc] = point_cointainment_cost(x,p_norm,I)
    n_I = length(I) ;
    x_vals = nan(n_I,1) ;
    
    for idx = 1:n_I
        x_vals(idx) = vecnorm(x(I{idx}),p_norm) ;
    end
    
    % compute cost
    [c,c_idx] = max(x_vals) ;
    
    % compute cost gradient
    gc = zeros(1,length(x)) ;
    gc(I{c_idx}) = p_norm*vecnorm(x(I{c_idx}),p_norm)*(x(I{c_idx}).^(p_norm-1)) ;
end