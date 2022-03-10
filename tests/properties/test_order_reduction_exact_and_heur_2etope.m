%% description
% This script tests order reduction for a 2-ellipsotope starting from some
% given number of generators and then cutting down to a user-specified
% number.
%
% See also: figure_order_reduction_2_etopes.m,
% test_order_reduction_high_dim.m
%
% Authors: Shreyas Kousik
% Created: 5 Mar 2022
% Updated: --
clear ; clc ; close all ;
%% user parameters
% rng seed
rng(1000)

% original etope properties
n_dim = 2 ; % leave this as 2 for plotting
n_gen = 20 ; % default is 20
n_con = 3 ; % default is 3
n_I = 3 ; % shrug

% reduction parameters
n_gen = 12 ; % default is 12

%% automated from here
% make random ellipsotope
p_norm = 2 ;
[E,c,G,A,b,I] = make_random_ellipsotope(p_norm,n_dim,n_gen,n_con,n_I) ;