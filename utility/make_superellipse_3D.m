function varargout = make_superellipse_3D(n,p,r,c)
% E = make_superellipse_3D(n,p,r,c)
% [F,V] = make_superellipse_3D(n,p,r,c)
%
% Make points on the boundary of the superellipse defined by the equation
%
%   (x_1 - c_1)^p + (x_2 - c_2)^p + (x_3 - c_3)^p = r^p
%
% where (x_1,x_2) are 2-D coordinates and c = (c_1,c_2) is a center point.
%
% All inputs are optional. The default input values are:
%   p = 2, r = 1, c = zeros(2,1), n = 100.
%
% The optional second output is a list of vertices so that you can plot a
% filled object with patch('faces',F,'vertices',V,other_patch_args)
%
% USAGE EXAMPLE;
%   [F,V] = make_superellipse_3D(4,2,zeros(3,1),200)
%   figure(1) ; axis equal ; hold on ; grid on ; view(3) ;
%   patch('faces',F,'vertices',V,'facealpha',0.1)
%
% See also: make_superellipse_2D
%
% Authors: Shreyas Kousik
% Created: 05 Mar 2021
% Updated: 17 Mar 2021

    if nargin < 1
        p = 2 ;
    end

    if nargin < 2
        r = 1 ;
    end
    
    if nargin < 3
        c = zeros(3,1) ;
    end

    if nargin < 4
        n = 100 ;
    end
    
    % create points using 2-norm
    E = make_sphere(1,zeros(3,1),n) ;
    
    % dilate the points using the given p-norm and shift the points to c
    N_E = vecnorm(E,p) ;
    E = r.*E./repmat(N_E,3,1) + repmat(c,1,size(E,2)) ;
    
    % create output
    if nargout == 1
        varargout = {E} ;
    elseif nargout == 2
        E_1 = E(1,:)' ;
        E_2 = E(2,:)' ;
        E_3 = E(3,:)' ;
        F = convhull(E_1,E_2,E_3) ;
        varargout = {F, E'} ;
    end
end