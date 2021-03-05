function [E,V] = make_superellipse_2D(r,p,n)
% E = make_superellipse_2D(r,p,n)
% [F,V] = make_superellipse_2D(r,p,n)
%
% Make points on the boundary of the superellipse defined by the equation
%
%   (x_1)^p + (x_2)^p + r^p = 1
%
% The optional second input is the number of points on the superellipse
% boundary (i.e., the length of the output E)
%
% The optional second output is a list of vertices so that you can plot a
% filled object with patch('faces',F,'vertices',V,other_patch_args)
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
        r = 1 ;
    end

    if nargin < 2
        p = zeros(2,1) ;
    end

    if nargin < 3
        n = 100 ;
    end
    
    
    if r < 0 || r > 1
        error('Please pick r in [0,1)')
    end
    
    if r == 1
        E = zeros(2,1) ;
    else
        % get max possible value of x_1 such that x_2 is real
        rp = r^p ;
        x_1_max = abs((1 - rp).^(1/p)) ; % setting x_2 = 0
        
        % make vector of points in first coord
        x_1_vec = linspace(-x_1_max,x_1_max,n/2) ;
        
        % make vector of points in second coord
        x_2_vec_up = (1 - r^p - (x_1_vec).^p).^(1/p) ;
        x_2_vec_dn = -x_2_vec_up ;
        E = [x_1_vec, x_1_vec(end-1:-1:1) ;
            x_2_vec_up, x_2_vec_dn(end-1:-1:1)] ;
    end


    if nargout > 1
        V = E' ;
        E = [1:n, 1] ;
    end
end