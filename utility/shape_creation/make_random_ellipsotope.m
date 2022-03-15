function [E,c,G,A,b,I] = make_random_ellipsotope(p_norm,n_dim,n_gen,n_con,n_I)
% E = make_random_ellipsotope(p_norm,n_dim,n_gen,n_con)
% [E,c,G,A,b,I] = make_random_ellipsotope(p_norm,n_dim,n_gen,n_con)
% make_random_ellipsotope(p_norm,n_dim,n_gen,n_con,I)
%
% Make a random ellipsotope. It does what it says on the box! The last
% input can also be an index set.
%
% Authors: Shreyas Kousik
% Created: 28 Apr 2021
% Updated: 15 Mar 2022 (added index set as optional input)

%% set default inputs
if nargin < 1
    p_norm = 2 ;
end

if nargin < 2
    n_dim = 2 ;
end

if nargin < 3
    n_gen = rand_int(2,10) ;
end

if nargin < 4
    n_con = rand_int(1,floor(n_gen/2)) ;
end

if nargin < 5
    n_I = rand_int(1,n_gen,n_gen/4,n_gen/4) ;
end

%% make random properties
c = 2*rand(n_dim,1) ;
G = 2*rand(n_dim,n_gen) - 1 ;
A = 2*rand(n_con,n_gen) - 1 ;
b = 0.5*rand(n_con,1) - 0.5 ;

if iscell(n_I)
    I = n_I ; 
    check_index_set_validity(I) ;
else
    I = make_random_index_set(n_gen,n_I) ;
end

%% create output
E = ellipsotope(p_norm,c,G,A,b,I) ;

end