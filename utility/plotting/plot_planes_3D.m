function h = plot_planes_3D(A,b,scale,varargin)
% h = plot_linear_subspace_3D(A,b)
% h = plot_linear_subspace_3D(A,b,scale)
% h = plot_linear_subspace_3D(A,b,scale,patch_args_in)
%
% Plot the planes defined by Ax = b.
%
% Authors: Shreyas Kousik
% Created: 9 Apr 2021
% Updated: nup

% get sizes of things
[n_planes,n_dim] = size(A) ;

% sanity check
if n_dim ~= 3 || length(b) ~= n_planes
    error(['Please make sure you pass in an n-by-3 array A and an n-by-1 ',...
        'array b.'])
end

% set the scale
if nargin < 3
    scale = 2.5 ;
end

% plot 'em and collect the handles
h = [] ;
for idx = 1:n_planes
    % get null space
    K = null(A(idx,:)) ;
    t = (A(idx,:))\b(idx) ;
    
    % plot plane
    F_H = [1 2 3 4 1] ;
    KK = [K' ; -K'] ;
    V_H = scale.*KK + repmat(t',size(KK,1),1) ;
    h_H = patch('faces',F_H,'vertices',V_H,...
        'facealpha',0.1','edgealpha',0,'facecolor','r',... % default arguments
        varargin{:}) ;
    
    % handle
    h = [h, h_H] ;
end

if nargout < 1
    clear h
end
end