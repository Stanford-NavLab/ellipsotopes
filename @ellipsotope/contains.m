function [out,value] = contains(E,p_test)
% [out,value] = E.contains(p)
%
% Check if the point p is inside the ellipsotope E, meaning the output
% value is <= 1.
%
% See also: ellipsotope.isempty
%
% Authors: Shreyas Kousik
% Created: 26 Apr 2021
% Updated: 27 Apr 2021

    % get etoproperties
    [p_norm,c,G,A,b,I] = E.get_properties() ;
    
    % compute point containtment constraint matrices
    A_eq = [G ; A] ;
    b_eq = [p_test - c ; b] ;

    % set up program cost
    cost = @(x) E.cost_for_emptiness_check(x,p_norm,I) ;
    x_0 = pinv(A_eq)*b_eq ;
    
    % set up options
    options = optimoptions('fmincon','Display','off',...
        'CheckGradients',false,... % useful to set to true sometimes...
        'SpecifyObjectiveGradient',true) ;
    
    % run optimization, woo!
    [~,value] = fmincon(cost,x_0,[],[],A_eq,b_eq,[],[],[],options) ;

    % output boolean
    out = value <= 1 ;
end