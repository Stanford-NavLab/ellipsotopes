function n = n_constraints(E)
% n = n_constraints(E)
%
% Return the number of constrats of the ellipsotope E
%
% Authors: Shreyas Kousik
% Created: 19 July 2021
% Updated: 20 July 2021 (fixed a big bug)
    n = size(E.constraint_A,1) ;
end