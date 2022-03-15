%% description
% This script demonstrates checking if a constrained zonotope is empty with
% an LP and a QP.
%
% Authors: Shreyas Kousik
% Created: 2 Jun 2021
% Updated: not yet
clear ; clc ;
%% user parameters
% rng seed
rng(0)

% zonotope spec
n_dim = 2 ;
n_gen = 10 ;
n_con = 3 ;

% whether or not the ellipsotope is empty
flag_empty = true ;

% run timing
flag_time_methods = true ;

%% automated from here
% make random ellipsotope (we'll just ignore the p_norm for the constrained
% zonotope that we care about)
[E,c,G,A,~,I] = make_random_ellipsotope(2,n_dim,n_gen,n_con) ;

% set b to be feasible or not
if flag_empty
    b = 10.*ones(n_con,1) ;
else
    b = zeros(n_con,1) ;
end

%% test LP
% cost
f_cost = [zeros(n_gen,1); 1] ;

% inequality cons
A_ineq = [-eye(n_gen), -ones(n_gen,1) ; eye(n_gen), -ones(n_gen,1)] ;
b_ineq = zeros(2*n_gen,1) ;

% equality cons
A_eq = [A, zeros(n_con,1)] ;
b_eq = b ;

% run LP
[~,value] = linprog(f_cost,A_ineq,b_ineq,A_eq,b_eq) ;
chk_empty_LP = value > 1

% time it
timeit(@() linprog(f_cost,A_ineq,b_ineq,A_eq,b_eq)) 

%% test QP
% cost
H = A'*A ;
f = (-2*b'*A)' ;
btb = b'*b ;

lb = -ones(n_gen,1) ;
ub = ones(n_gen,1) ;

% run QP
[~,value] = quadprog(H,f,[],[],[],[],lb,ub) ;
value = value + btb ;

chk_empty_QP = abs(value) > 1e-10 

% time it
timeit(@() quadprog(H,f,[],[],[],[],lb,ub)) 