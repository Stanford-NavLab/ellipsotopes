function C = make_unit_superellipse_ND(p,d,n)
% C = make_unit_superellipse_ND(p,d,n)
%
% Given a d-dimensional unit superellipse defined by the p-norm, return n
% points on the superellipse.
%
% All of the inputs are optional, with defaults:
%   p = 2, d = 2, n = 100
%
% The output is a d-by-n array of random points on the superellipse,
% sampled roughly uniformly.
%
% Authors: Shreyas Kousik
% Created: 15 Mar 2021

    %% set default inputs
    if nargin < 1
        p = 2 ;
    end

    if nargin < 2
        d = 2 ;
    end

    if nargin < 3
        n = 100 ;
    end
    
    n_2 = 2^d ;
    if n_2 < n
        % generate points on n-D grid "corners" of hypercube if we can do
        % so without things blowing up
        P = make_grid(repmat([-1,1],1,d),repmat(2,1,d)) ;
        
        % add additional points to fill things in
        P = [P, 2*rand(2,n_2 - n) - 1] ;
    else
        % generate n points
        P = 2.*rand(d,n) - 1 ;
    end
    
    % map all points to the p-sphere
    C = P./repmat(vecnorm(P,p),d,1) ;
end