function [E,V] = make_unit_superellipse_3D(p,n)


if nargin < 1
    p = zeros(2,1) ;
end

if nargin < 2
    n = 30 ;
end

% make points in all three coords
n = ceil(n/3) ;
X = make_grid([-1,1,-1,1,-1,1],[n n n]) ;

% evaluate which points obey the norm constraint
p_log = sum(X.^p,1) <= 1 ;

% get boundary of points
X = X(:,p_log) ;
K = boundary(X') ;
X = X(:,K) ;

if nargout == 1
    E = X ;
else
    V = X' ;
    E = boundary(V) ;
end
end