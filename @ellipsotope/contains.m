function [result,value] = contains(E,p_test)
% [result,value] = E.contains(p)
%
% Check if the point p is inside the ellipsotope E, meaning the output
% value is <= 1.
%
% See also: ellipsotope.isempty
%
% Authors: Shreyas Kousik
% Created: 26 Apr 2021
% Updated: 01 Jun 2021 (used intersection and E.isempty())

    % get p_norm
    p_norm = E.get_properties() ;
    
    % create an ellipsotope with just the center point
    E_p = ellipsotope(p_norm,p_test,[]) ;
    
    % compute intersection
    E_int = E & E_p ;
    
    % compute output
    if nargout > 1
        [out,value] = E_int.isempty(true) ;
    else
        out = E_int.isempty(false) ;
    end
    
    result = ~out ;
end