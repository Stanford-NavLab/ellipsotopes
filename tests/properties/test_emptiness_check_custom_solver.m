%% description
% This script builds a solver for checking if an ellipsotope is empty.
%
% Authors: Shreyas Kousik
% Created: 27 Apr 2021
% Updated: 28 Apr 2021
clear ; clc ;

%% user parameters
% rng seed
rng(0)

% ellipsotope
p_norm = 2 ;
c = zeros(2,1) ;
G = 2*rand(2,4) - 1 ;
A = rand(1,4) ;
b = 10 ; % nonempty
% b = 100 ; % empty
I = {1:2,3,4} ;

% solver timeout
t_max = 1 ; % seconds

%% automated from here
% sizes of things
n_I = length(I) ;
n_x = size(G,2) ;

% get nullspace of constraint
K = null(A) ;
K = K ./ vecnorm(K,2) ; % make sure the nullspace vectors are unit length

% create initial feasible point on linear constraint
t = pinv(A)*b ;

% preallocate values to max over
c_vals = nan(n_I,1) ;

%% run solver
% this flag stays truuuuu
flag_iter = true ;

% current iterate
x = t ;

start_tic = tic ;

while flag_iter && (toc(start_tic) < t_max)
    % compute the current cost
    for idx = 1:n_I
        c_vals(idx) = sum(x(I{idx}).^p_norm) ;
    end
    
    % compute cost
    [c,c_idx] = max(c_vals) ;
    
    if c <= 1
        flag_iter = true ;
        break
    else
        J_max = I{c_idx} ;

        % compute cost (sub)gradient
        c_grad = zeros(n_x,1) ;
        c_grad(J_max) = p_norm*(x(J_max).^(p_norm-1)) ;

        % compute cost (sub)Hessian
        c_hess = zeros(n_x) ;
        J_max_ind = sub2ind([n_x,n_x],J_max,J_max) ;
        c_hess(J_max_ind) = (p_norm*(p_norm-1))*(x(J_max).^(p_norm-2)) ;
        
        % compute descent direction
        dx = pinv(c_hess)*c_grad ;

        % project descent direction onto linear subspace
        dx_dots = K'*(dx - t) ;
        dx = K*dx_dots ;
        
        % take a step
        x = x - dx ;
    end
end

%% compare against implemented class method