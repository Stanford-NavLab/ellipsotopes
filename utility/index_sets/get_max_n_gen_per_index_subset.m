function n_gen_max = get_max_n_gen_per_index_subset(I)
% n_gen_max = get_max_n_gen_per_index_subset(I)
%
% Get the maximum number of generators indexed together by any of the index
% subsets J within the index set I
%
% Authors: Shreyas Kousik
% Created: 14 March 2022
% Updated: nah
    n_per_J = cellfun(@(J) length(J),I) ;
    n_gen_max = max(n_per_J) ;
end