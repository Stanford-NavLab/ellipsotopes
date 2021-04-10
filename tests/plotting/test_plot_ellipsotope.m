%% description
% This script tests plotting a (generalized, constrained) ellipsotope
% (without using the superclass!)
%
% TO DO:
%   - put plotting code into ellipsotope.m
%   - validate plotting code using the point containment lemma
%
% See also: test_projection_to_constraint_ball_product_boundary.m
%
% Authors: Shreyas Kousik
% Created: 9 Apr 2021
% Updated: nope
clear ; clc ;
%% user parameters
% rng seed
rng(0)

% ellipsotope definition (make is 2-D please)
p_norm = 2 ;
c = zeros(2,1) ;
G = 2*rand(2,4) - 1;
A = [-1 1 -1 1] ;
b = 0.5 ;
I = {[1,2],[3,4]} ;

%% automated from here
% TO DO: functionize these sanity checks within the ellipsotope.m class as
% a method

% sanity check the index set
[I_chk,n_I,n_dim] = check_index_set_validity(I) ;
if ~I_chk
    error(['The index set is not valid! ',...
        'It should contain each dimension of the hyperball product ',...
        'exactly once.'])
end

% get number of constraints
[n_con,n_dim_con] = size(A) ;

% sanity check the constraints
if n_dim_con ~= n_dim
    error(['The constraints defined by (A,b) must have as many columns ',...
        'as there are dimensions in the index set (i.e., the hyperball ',...
        'space of the ellipsotope.'])
end

%% plotting setup
tic
% create a bunch of points in the ball space
n_P = 1000 ;
P = 2*rand(n_dim,n_P) - 1 ;

% project points to intersection of ball product and linear subspace
[P,n_P] = project_points_to_ball_product_and_linear_subspace(P,p_norm,A,b,I) ;


% map points to workspace
P = repmat(c,1,n_P) + G*P ;
K = convhull(P') ;
P = P(:,K) ;
P = points_to_CCW(P) ;
n_P = size(P,2) ;
toc

%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ;

% plot ellipsotope
patch('faces',[1:n_P,1],'vertices',P','facecolor','b','facealpha',0.1,...
    'linewidth',1.5,'edgecolor','b')