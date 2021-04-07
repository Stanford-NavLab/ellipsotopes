%% description
% This script tests doing a line search in a random direction to find the
% point on the line where the generalized norm is equal to 1. The purpose
% of this is to find points on the boundary of a product of balls in any
% given direction, to enable plotting generalized ellipsotopes.
%
% Authors: Shreyas Kousik
% Created: 6 Apr 2021
% Updated: nup
clear;clc
%% user parameters
% rng seed
rng(0) ;

% p norm
p_norm = 2 ;

% index set
I = {1,[2,3]} ;

%% automated from here

%% plotting setup

%% plotting