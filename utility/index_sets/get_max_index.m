function v = get_max_index(I)
% v = get_max_index(I)
%
% For an index set I, get the maximum index (should be the same as the
% number of generators)
%
% Authors: Shreyas Kousik
% Created: 3 Mar 2022
% Updated: nope
    if isempty(I)
        v = 0 ;
    else
        v = max(cell2mat(I)) ;
    end
end