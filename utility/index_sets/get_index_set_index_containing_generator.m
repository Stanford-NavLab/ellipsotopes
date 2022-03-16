function idx_J = get_index_set_index_containing_generator(I,idx_gen)
% idx_J = get_index_set_index_containing_generator(I,idx_gen)
%
% For an index set I, return the index of J \in I for which J contains the
% integer idx_gen; that is idx_gen \in J = I{idx_J}
%
% Authors: Shreyas Kousik
% Created: 16 Mar 2022
% Updated: nup

    % get the index subsets containing idx_gen
    idx_log = cellfun(@(J) any(J == idx_gen),I) ;
    idxs_I = 1:length(I) ;
    idx_J = idxs_I(idx_log) ;
end