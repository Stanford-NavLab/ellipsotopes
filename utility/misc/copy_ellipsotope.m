function E_copy = copy_ellipsotope(E)
% E_copy = copy_ellipsotope(E)
%
% Does what it says on the lid.
%
% Authors: Shreyas Kousik
% Created: 14 Mar 2022
% Updated: nope
    [p,c,G,A,b,I] = get_properties(E) ;
    E_copy = ellipsotope(p,c,G,A,b,I) ;
end