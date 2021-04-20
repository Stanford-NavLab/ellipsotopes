%% description
% This script tests swapping the min/max of the point containment lemma; it
% turns out this doesn't work because Shreyas applied Von Neumann's minimax
% theorem incorrectly :)
%
% Authors: Shreyas Kousik
% Created: 20 Apr 2021
% Updated:
clear ; clc
%% user parameters
% rng seed
rng(0)

% point to test
% p_test = [1.0 ; 1.0] ; % this point is inside
p_test = [1.45; 0.0] ; % this point is barely outside

% ellipsotope
p_norm = 2 ;
c = zeros(2,1) ;
G = 2*rand(2,4) - 1 ;
A = [] ;
b = [] ;
I = {1:2,3,4} ;

% whether or not to time it
flag_time_check = false ;

%% automated from here
% get sizes of things
n_gen = size(G,2) ;
n_con = size(A,1) ;
n_I = length(I) ;

% compute point containtment constraint matrices
A_eq = [G ; A] ;
b_eq = [p_test - c ; b] ;

% optimization options
options = optimoptions('fmincon','Display','off','SpecifyObjectiveGradient',true) ;

%% min of maxs
% set up program cost
cost = @(x) point_cointainment_cost_minmax(x,p_norm,I) ;

% take initial guess
x_0 = pinv(A_eq)*b_eq ;

% solve problem
[x_opt_minmax,f_val_minmax] = fmincon(cost,x_0,[],[],A_eq,b_eq,[],[],[],options) ;

disp('minmax optimal cost:')
disp(f_val_minmax)

%% max of mins
% initial guess
x_0 = pinv(A_eq)*b_eq ;

% store solutions
f_val_maxmin = nan(1,n_I) ;

% iterate through index set and solve subproblems
for idx = 1:n_I
    % get current index set
    J = I{idx} ;
    
    % cost function
    cost = @(x) point_cointainment_cost_maxmin(x,p_norm,J) ;
    
    % solve subproblem
    [~,f_val_idx] = fmincon(cost,x_0,[],[],A_eq,b_eq,[],[],[],options) ;

    % save solution
    f_val_maxmin(idx) = f_val_idx ;
end

disp('maxmin optimal costs:')
disp(f_val_maxmin)

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
function [c,gc] = point_cointainment_cost_minmax(x,p_norm,I)
    n_I = length(I) ;
    x_vals = nan(n_I,1) ;

    for idx = 1:n_I
        x_vals(idx) = sum(x(I{idx}).^p_norm) ;
    end

    % compute cost
    [c,c_idx] = max(x_vals) ;

    % compute cost (sub)gradient
    gc = zeros(1,length(x)) ;
    gc(I{c_idx}) = p_norm*(x(I{c_idx}).^(p_norm-1)) ;
end

function [c,gc] = point_cointainment_cost_maxmin(x,p_norm,J)
    c = sum(x(J).^p_norm) ;
    gc = zeros(1,length(x)) ;
    gc(J) = p_norm.*(x(J)').^(p_norm-1) ;
end