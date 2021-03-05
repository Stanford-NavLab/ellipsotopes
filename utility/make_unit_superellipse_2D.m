function [E,V] = make_unit_superellipse_2D(p,n)
% E = make_unit_superellipse_2D(p,n)
% [F,V] = make_unit_superellipse_2D(p,n)
%
% Make points on the boundary of the superellipse defined by the equation
%
%   (x_1)^p + (x_2)^p = 1
%
% The optional second input is the number of points on the superellipse
% boundary (i.e., the length of the output E)
%
% The optional second output is a list of vertices so that you can plot a
% filled object with patch('faces',F,'vertices',V,other_patch_args)
%
% This uses the parametric Lam'e curve representation from here:
%   https://mathworld.wolfram.com/Superellipse.html
%
% USAGE EXAMPLE;
%   E = make_unit_superellipse_2D(4,200)
%   figure(1) ; axis equal ; hold on ; grid on ;
%   plot_path(E,'b-')
%
% Authors: Shreyas Kousik
% Created: 5 Mar 2021
% Updated: nah

    if nargin < 1
        p = zeros(2,1) ;
    end

    if nargin < 2
        n = 100 ;
    end
    
    % make vector of angles for positive quadrant
    n = ceil(n/4) ;
    t = linspace(0,pi/2,n) ;

    % make points for the positive quadrant
    X = [cos(t).^(2/p) ;
         sin(t).^(2/p) ] ;
     
    % make rotation matrices
    R = rotation_matrix_2D(pi/2) ;
    
    % create points
    E = [X, R*X(:,2:end), R^2*X(:,2:end), R^3*X(:,2:end)] ;

    if nargout > 1
        V = E' ;
        E = [1:n, 1] ;
    end
end