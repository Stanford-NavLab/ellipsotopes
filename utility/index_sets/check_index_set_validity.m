function [out,n_I,n_gen] = check_index_set_validity(I,G)
% out = check_index_set_validity(I)
% [out,n_I,n_gen] = check_index_set_validity(I,G)
%
% Given an ellipsotope index set, return True if it is valid and False if
% it is not. The generator matrix can be passed in as an optional argument.
% The optional outputs are the number of index subsets and the number of
% generators of the ellipsotope.
%
% Authors: Shreyas Kousik
% Created: 4 Apr 2021
% Updated: 15 Mar 2022 (fixed bug when G is passed in)

    % optimism
    out = true ;

    % number of index subsets
    n_I = length(I) ;
    
    % indices all in one vector
    D = cell2mat(I) ;

    % get the dimension of the ball product
    if nargin > 1
        n_gen = size(G,2) ;
    else
        n_gen = max(D) ;
    end
    
    % check!
    if (length(unique(D)) < length(D)) || ... % no repeated indices
            (length(D) ~= n_gen) || ... % correct length
            (max(D) > n_gen) || (min(D) ~= 1) % not shifted to wrong index
        out = false ;
    end
end
