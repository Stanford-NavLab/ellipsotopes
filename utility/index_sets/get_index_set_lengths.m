function L = get_index_set_lengths(I)
% L = get_index_set_lengths(I)
%
% For an index set of length n_I, get the length of each J \in I
%
% Authors: Shreyas Kousik
% Created: 16 Mar 2022
% Updated: not yet
    L = cellfun(@(J) length(J), I) ;
end