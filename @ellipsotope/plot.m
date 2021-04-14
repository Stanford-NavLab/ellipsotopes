function plot(E,varargin)
% plot(E)
% plot(E,'projdims',[dim1 dim2], other_input_args...)
% plot(E,'facecolor',color,'edgecolor',color,'facealpha',...)
%
% Plot the ellipsotope if it is 2-D. This creates a patch
% object, and updates the E.plot_handle property.
%
% To plot an ellipsotope, we do the following:
%   1. generate a bunch of points in the coefficient space
%   2. push those points to the boundary of the feasible coefficient space
%   3. send the points through the affine map defined by the ellipsotope's
%      center and generators
%   4. plot the convex hull of the mapped points
%
% See also: test_plot_ellipsotope.m
%
% Authors: Adam Dai and Shreyas Kousik
% Created: in days of yore
% Updated: 13 Apr 2021

%% prep/sanity check
% check the dimension
if E.dimension > 2
    warning(['Plotting not supported for > 2-D ellipsotopes!',...
        'Plotting a 2-D projection instead'])
    % TODO: select projection dims for user?
end

% check for projection dimensions
if nargin > 1 && strcmpi(varargin{1},'projdims')
    proj_dims = varargin{2} ;
    
    % clean up varargin
    if length(varargin) > 2
        varargin = varargin(3,:) ;
    end
end

% get important properties
p = E.p_norm ;
c = E.center ;
G = E.generators ;
I = E.index_set ;
d = E.dimension ;
d_B = size(G,2) ; % dimension of coefficient space

% set proj dims if dimension is greater than 2
if ~exist('proj_dims','var')
    if d <= 2
        proj_dims = 1:d ;
    else
        proj_dims = 1:2 ;
        warning(['Only plotting first ',num2str(d),' dimensions!'])
    end
end
G = G(proj_dims,:) ;

% check if E is basic, which allows us to reduce the
% generator matrix nicely
if E.is_basic() && (p == 2)
    % check if E is reduced
    if ~E.is_reduced()
        G = reduce_ellipsotope_generator_matrix(G) ;
    end
end

% set the index set if it is empty
if isempty(I)
    I = {1:d_B} ; % constrain ALL the coefficients!
end

%% plotting setup
% STEP 1: generate points in coefficient space
% if the ellipsotope is basic...
    % switch how we generate points based on the coefficient space
    % dimension (using make_superllipse_2D, _3D, or _ND)
    
% else if...
    % generate a bunch of random points 
    
    % project points to ball and linear subspace boundary (the existing
    % function will handle index sets and empty linear subspaces properly,
    % it turns out)
    
% STEP 2: map points to ellipsotope workspace

% STEP 3: take convex hull of points in workspace

%% plotting code here
% now we actually call patch, hehe

%% OLD CODE FROM HERE ON
% generate a bunch of points in the unit hypercube space and
% plot the resulting object; we use "B" for "beta" which
% we've used in our paper's notation for the ellipsotope
% coefficients
switch d_B
    case 2
        n_plot = 100 ;
        B = make_superellipse_2D(p,1,zeros(2,1),n_plot) ;
    case 3
        n_plot = 1000 ;
        B = make_superellipse_3D(p,1,zeros(3,1),n_plot) ;
    otherwise
        % get points on superellipse
        n_plot = 10000 ; % we can be smarter than this I guess
        B = make_unit_superellipse_ND(p,d_B,n_plot) ;
end

% check if we need to project the points onto a constraint
if E.is_constrained()
    % get the constraints
    A = E.constraint_A ;
    b = E.constraint_b ;
    
    % get a basis for the constraint space
    N = null(A) ;
    b_offset = pinv(A)*b ;
    
    % project points onto the plane
    B = pinv(N*N')*(B - b_offset) + b_offset ;
    
    % add some more points in the nullspace (shrugs)
    n_extra = 10000 ;
    B_extra = 4*rand(size(N,2),n_extra) - 2 ;
    B = [B, N*B_extra + b_offset] ;
    
    % keep the points that obey the norm
    for idx = 1:length(E.index_set)
        N_log = vecnorm(B(E.index_set{idx},:),p) <= 1 ;
        B = B(:,N_log) ;
    end
end

% map points to the ellipsotope's proj dims (recall that the
% generators have already been projected)
P = G*B + repmat(c,1,size(B,2)) ;

% get the convex hull of the points
ch = convhull(P') ;

% set default plot inputs
patch_options = [{'facecolor','b','edgecolor','b',...
    'facealpha',0.1,'edgealpha',1}, varargin] ;

% plot!
switch length(proj_dims)
    case 1
        error('1-D plot is not implemented yet!')
    case 2
        P = P(:,ch) ;
        P = points_to_CCW(P) ;
        F = [1:size(P,2),1] ;
        h = patch('faces',F,'vertices',P',patch_options{:}) ;
    case 3
        h = trisurf(ch,P(1,:)',P(2,:)',P(3,:)',patch_options{:}) ;
end

% finish up
E.plot_handle = h ;
end