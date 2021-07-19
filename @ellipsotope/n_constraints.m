function n = n_constraints(E)
% n = n_constraints(E)
%
% Return the number of constrats of the ellipsotope E
%
% Authors: Shreyas Kousik
% Created: 19 July 2021
% Updated: --
    n = size(E.constraint_A,2) ;
end