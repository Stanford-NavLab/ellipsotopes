%% description
% This script demonstrates order reduction for a 2-etope where we get rid
% of constraints by leveraging the fact that the intersection of a
% hyperplane with a 2-norm ball is an affine map of a lower-dimensional
% 2-norm ball
%
% See also: test_ball_slice_lemma.m
%
% Authors: Shreyas Kousik
% Created: 13 Jul 2021
% Updated: -
clear ; clc

%% user parameters
% rng seed
rng(0) ;

% original tope
n_gen = 10 ;
n_con = 3 ;

%% automated from here
% make original etope
[E,c,G,A,b,I] = make_random_ellipsotope(2,2,n_gen,n_con,1) ;

% create new tope
E_rdc = reduce_constrained_2_etope(E) ;

%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ;
plot(E)
plot(E_rdc,'color','r','linestyle','--','linewidth',3)

legend('orig','reduced')

set_plot_fontsize(15)