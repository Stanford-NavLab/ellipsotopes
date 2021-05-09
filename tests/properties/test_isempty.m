%% description
% This script tests the "isempty" emptiness check function for ellipsotopes
%
% Authors: Adam Dai
% Created: 9 May 2021
% Updated: 
clear ; clc ; 
%% user parameters
% rng seed
rng(0)

% ellipsotope
p_norm = 2 ;
n_dim = 2 ;
n_gen = 10 ;
n_con = 1 ;


%% automated from here
% make random etope
E = make_random_ellipsotope(p_norm,n_dim,n_gen,n_con) ;

% emptiness check 
[out,value] = isempty(E);
