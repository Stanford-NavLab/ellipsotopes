%% description
% This script compares the solve speed of different methods for checking if
% an ellipsotope is empty. It looks like method 2 is quite a bit faster.
%
% Authors: Shreyas Kousik
% Created: 28 Apr 2021
% Updated: 1 Jun 2021 (fixed gradient for method 2)
clear ; clc ; 
%% user parameters
% rng seed
rng(0)

% ellipsotope
p_norm = 2 ;
n_dim = 5 ;
n_gen = 10 ;
n_con = 2 ;
n_idx = 4 ;

% whether or not the ellipsotope is empty
flag_empty = false ;

% run timing
flag_time_methods = true ;

%% automated from here
% make random etope
[E,c,G,A,~,I] = make_random_ellipsotope(p_norm,n_dim,n_gen,n_con,n_idx) ;

% set b to be feasible or not
if flag_empty
    b = 10.*ones(n_con,1) ;
else
    b = 0.25*rand(n_con,1) - 0.125 ;
end

E.constraint_b = b ;
 
%% problem setup
% initial guess
x_0 = pinv(A)*b ;

% fmincon options
options = optimoptions('fmincon','Display','off',...
    'SpecifyObjectiveGradient',true,...
    'SpecifyConstraintGradient',true) ;

%% emptiness check v1
cost = @(x) cost_for_emptiness_check_v1(x,p_norm,I) ;
[x_opt_1,f_val_1] = fmincon(cost,x_0,[],[],A,b,[],[],[],options) ;

if f_val_1 < 1
    disp('Ellipsotope is non-empty by method 1!')
else
    disp('Ellipsotope is empty by method 1!')
end

if flag_time_methods
    t_avg_1 = timeit(@() fmincon(cost,x_0,[],[],A,b,[],[],[],options)) ;
    disp(['Method 1 takes ',num2str(1000*t_avg_1,'%0.2f'),' ms to solve on average'])
end

%% emptiness check v2
cost_2 = @(x) cost_for_emptiness_check_v2(x,A,b) ;
nonlcon = @(x)  nonlcon_for_emptiness_check(x,p_norm,I) ;
[x_opt_2,f_val_2] = fmincon(cost_2,x_0,[],[],A,b,[],[],nonlcon,options) ;

if abs(f_val_2) < 1e-10
    disp('Ellipsotope is non-empty by method 2!')
else
    disp('Ellipsotope is empty by method 2!')
end

if flag_time_methods
    t_avg_2 = timeit(@() fmincon(cost_2,x_0,[],[],[],[],[],[],nonlcon,options)) ;
    disp(['Method 2 takes ',num2str(1000*t_avg_2,'%0.2f'),' ms to solve on average'])
end

%% helper functions
function [c,gc] = cost_for_emptiness_check_v1(x,p_norm,I)
    % get number of index subsets
    n_I = length(I) ;
    
    % preallocate values to max over
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

function [c,gc] = cost_for_emptiness_check_v2(x,A,b)
    c = sum((A*x - b) .^2) ;
    gc = (2*(A'*A)*x)' - 2*b'*A ;
end

function [c,ceq,gc,gceq] = nonlcon_for_emptiness_check(x,p_norm,I)
    % get sizes of things
    n_x = length(x) ;
    n_I = length(I) ;
    
    % preallocate constraint and gradient
    c = zeros(n_I,1) ;
    gc = zeros(n_x,n_I) ;
    
    % compute constraint and gradient
    for idx = 1:n_I
        J = I{idx} ;
        c(idx) = sum(x(J).^p_norm) - 1 ;
        gc(J,idx) = p_norm*x(J).^(p_norm-1) ;
    end
    
    % equality constraint outputs
    ceq = [] ;
    gceq = [] ;
end