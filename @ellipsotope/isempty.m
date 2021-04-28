function [out,value] = isempty(E)
% [out,value] = isempty(E)
%
% This function checks if the ellipsotope E is empty by solving a convex
% program. The output is TRUE if the ellipsotope is EMPTY.
%
% See also: ellipsotope.contains
%
% Authors: Shreyas Kousik
% Created: 27 Apr 2021
% Updated: nuuu

    % get etoproperties
    [p_norm,~,~,A,b,I] = E.get_properties() ;

    if isempty(A)
        % no constraints = not empty
        out = false ;
        value = 0 ;
    else
        % set up program cost
        cost = @(x) E.cost_for_emptiness_check(x,p_norm,I) ;
        x_0 = pinv(A)*b ;
        
        % set up options
        options = optimoptions('fmincon','Display','off',...
            'CheckGradients',false,... % useful to set to true sometimes...
            'SpecifyObjectiveGradient',true) ;
        
        % run optimization, woo!
        [~,value] = fmincon(cost,x_0,[],[],A,b,[],[],[],options) ;
        
        % output boolean
        out = value > 1 ;
    end
end