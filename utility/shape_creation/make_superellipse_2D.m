function [out_1,out_2] = make_superellipse_2D(varargin)
% E = make_superellipse_2D(n,p,r,c)
% [F,V] = make_superellipse_2D(n,p,r,c)
%
% Make points on the boundary of the superellipse defined by the equation
%
%   (x_1 - c_1)^p + (x_2 - c_2)^p = r^p
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
%   E = make_superellipse_2D(4,2,rand(2,1),200)
%   figure(1) ; axis equal ; hold on ; grid on ;
%   plot_path(E,'b-')
%
% See also: make_superellipse_3D
%
% Authors: Shreyas Kousik
% Created: 05 Mar 2021
% Updated: 17 Mar 2021
    
    % defaults 
    n = 100; p = 2; r = 1; c = zeros(2,1);
    
    if nargin > 0
        n = varargin{1} ;
    end

    if nargin > 1
        p = varargin{2} ;
    end
    
    if nargin > 2
        r = varargin{3} ;
    end

    if nargin > 3
        c = varargin{4} ;
    end
    
    % create points using 2-norm
    E = make_circle(1,n) ;
    
    % dilate the points using the given p-norm and shift the points to c
    N_E = vecnorm(E,p) ;
    E = r.*E./N_E + repmat(c,1,n);
    
    % create output
    if nargout == 1
        out_1 = E ;
    elseif nargout == 2
        out_1 = [1:n,1] ;
        out_2 = E' ;
    end
end