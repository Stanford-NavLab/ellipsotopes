%% description
% This script identifies component ellipsotopes by checking for block
% diagonal entries in the constraint matrix according to the index set.
%
% Authors: Shreyas Kousik
% Created: 13 July 2021
% Updated: not yet
clear ; clc ;

%% automated from here
rng(0)
E_1 = make_random_ellipsotope(2,2,4,1,1) ;
E_2 = make_random_ellipsotope(2,2,4,1,1) ;
E = E_1 + E_2 ;

% get properties
[p,c,G,A,b,I,n_dim,n_gen,n_con,n_I] = get_properties(E)

% initialize output
idxs = [] ;

% iterate through index set and check for component etopes
for idx = 1:n_I
    % get current indices
    J = I{idx} ;
    
    % get all constraints with nonzero values at the current indices
    A_row_log = any(A(:,J) ~= 0,2) ;
    
    % check if all other entires in the current constraints are zeros
    A_temp = A ;
    A_temp(:,J) = [] ;
    if all(all(A_temp(A_row_log,:),2) == 0)
        idxs = [idxs, idx] ;
    end
end

%%
figure(1) ; clf ; axis equal ; hold on ; 

plot(E)