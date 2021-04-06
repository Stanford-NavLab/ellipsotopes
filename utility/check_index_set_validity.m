function [out,n_I,n_dim] = check_index_set_validity(I)
% out = check_index_set_validity(I)
%
% Given an ellipsotope index set, return True if it is valid and False if
% it is not.
%
% Authors: Shreyas Kousik
% Created: 4 Apr 2021
% Updated: nup

    % optimism
    out = true ;

    % number of index subsets
    n_I = length(I) ;
    
    % indices all in one vector
    D = cell2mat(I) ;

    % get the dimension
    n_dim = max(D) ;

    if length(unique(D)) < length(D) || length(D) < n_dim
        out = false ;
    end
end
