%% description
% This script compares the solve speed of different methods for checking if
% an ellipsotope is empty. It looks like method 1 is faster.
%
% Authors: Shreyas Kousik
% Created: 28 Apr 2021
% Updated: nope
clear ; clc ; 
%% user parameters
% rng seed
rng(0)

% ellipsotope
p_norm = 2 ;
n_dim = 2 ;
n_gen = 10 ;
n_con = 1 ;

% whether or not the ellipsotope is empty
flag_empty = true ;

%% automated from here
% make random etope
E = make_random_ellipsotope(p_norm,n_dim,n_gen,n_con) ;
[p_norm,c,G,A,~,I] = get_properties(E) ;

% set b to be feasible or not
if flag_empty
    b = 100.*ones(n_con,1) ;
else
    b = zeros(n_con,1) ;
end

%% problem setup
% initial guess
x_0 = pinv(A)*b ;

% fmincon options
options = optimoptions('fmincon','Display','off',...
    'SpecifyObjectiveGradient',true,...
    'SpecifyConstraintGradient',true) ;

%% emptiness check v1
cost = @(x) cost_for_emptiness_check(x,p_norm,I) ;
[x_opt_1,f_val_1] = fmincon(cost,x_0,[],[],A,b,[],[],[],options) ;

% timeit(@() fmincon(cost,x_0,[],[],A,b,[],[],[],options))

%% emptiness check v2
cost_2 = @(x) cost_for_emptiness_check_v2(x) ;
nonlcon = @(x)  nonlcon_for_emptiness_check(x,p_norm,I) ;
[x_opt,f_val] = fmincon(cost_2,x_0,[],[],A,b,[],[],nonlcon,options) ;

% timeit(@() fmincon(cost_2,x_0,[],[],A,b,[],[],nonlcon,options))

%% helper functions
function [c,gc] = cost_for_emptiness_check(x,p_norm,I)
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

function [c,gc] = cost_for_emptiness_check_v2(x)
    c = sum(x.^2) ;
    gc = 2*x(:) ;
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