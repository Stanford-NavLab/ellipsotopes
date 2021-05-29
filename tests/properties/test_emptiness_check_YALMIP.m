%% description
% This script constructs an ellipsotope emptiness check as a feasibility
% problem using YALMIP.
%
% Unfortunately, this runs slower than fmincon for the always-feasible
% problem.
%
% See also: test_emptiness_check.m
%
% Authors: Shreyas Kousik
% Created: 28 May 2021
% Updated: not yet
clear ; clc ;
%% user parameters
% rng seed
rng(0)

% ellipsotope
p_norm = 4 ;
n_dim = 3 ;
n_gen = 15 ;
n_con = 1 ;

% whether or not the ellipsotope is empty
flag_empty = false ;

%% automated from here
disp('Creating random ellipsotope') 

% make random etope
E = make_random_ellipsotope(p_norm,n_dim,n_gen,n_con) ;
[p_norm,c,G,A,~,I] = get_properties(E) ;

% set b to be feasible or not
if flag_empty
    b = 100.*ones(n_con,1) ;
else
    b = zeros(n_con,1) ;
end

%% create YALMIP program
disp('Creating YALMIP problem')

% create decision variable
x = sdpvar(n_gen,1) ;

% create linear constraint
cons = A*x == b ;

% create norm constraints
n_I = length(I) ;
for idx = 1:n_I
    J = I{idx} ;
    cons = [cons, sum(x(J).^p_norm) - 1] ;
end

%% solve YALMIP problem
disp('Solving YALMIP problem')

tic
optimize(cons) ;
toc

%% solve with usual emptiness check
disp('Solving with regular emptiness check')

tic
E.isempty() ;
toc