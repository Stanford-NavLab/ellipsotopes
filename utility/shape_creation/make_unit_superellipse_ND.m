function C = make_unit_superellipse_ND(varargin)
% C = make_unit_superellipse_ND(n,p,d)
%
% Given a d-dimensional unit superellipse defined by the p-norm, return n
% points on the superellipse.
%
% All of the inputs are optional, with defaults:
%   n = 100, p = 2, d = 2
%
% The output is a d-by-n array of random points on the superellipse,
% sampled roughly uniformly.
%
% Authors: Shreyas Kousik
% Created: 15 Mar 2021
% Updated: 17 Mar 2021

    %% set default inputs
    % defaults 
    n = 100; p = 2; d = 2;

    if nargin > 0
        n = varargin{1} ;
    end

    if nargin > 1
        p = varargin{2} ;
    end
    
    if nargin > 2
        d = varargin{3} ;
    end
    
    n_2 = 2^d ;
    if n_2 < n
        % generate points on n-D grid "corners" of hypercube if we can do
        % so without things blowing up
        P = make_grid(repmat([-1,1],1,d),repmat(2,1,d)) ;
        
        % add additional points to fill things in (hopefully)
        P = [P, 2*rand(d,n - n_2) - 1] ;
    else
        % generate n points
        P = 2.*rand(d,n) - 1 ;
    end
    
    % map all points to the p-sphere
    C = P./repmat(vecnorm(P,p),d,1) ;
end