function result = subset(E1,E2)
% result = subset(E1,E2)
% result = E1.subset(E2)
%
% Check if ellipsotope E1 is a subset of ellipsotope E2 (i.e. if E2
% contains E1)
%
% Authors: Adam Dai
% Created: 30 Jun 2021
% Updated: 

    % compute intersection
    E_int = E1 & E2 ;
    
    % check if intersection equals E1
    
    
    result = ~out ;
end