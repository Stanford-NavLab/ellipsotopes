function I_out = shift_index_set(I_in,n_shift)
% I_out = shift_index_set(I_in,n_shift)
%
% For each J \in I, output an index set containing each J + n_shift
%
% Authors: Shreyas Kousik
% Created: 21 Apr 2021
% Updated: nah
    I_out = cellfun(@(x) x + n_shift,I_in,'UniformOutput',false) ;
end