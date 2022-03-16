function E_lift = lift(E)
% E_lift = lift(E)
%
% Lift E to be an indexed-but-not-constrained ellipsotope.
%
% Authors: Shreyas Kousik
% Created: 15 Mar 2022
% Updated: nah

    [p,c,G,A,b,I] = get_properties(E) ;
    
    % lift
    G_lift = [G ; A] ;
    c_lift = [c ; -b] ;
    
    % make new outputope
    E_lift = ellipsotope(p,c_lift,G_lift,[],[],I) ;
end