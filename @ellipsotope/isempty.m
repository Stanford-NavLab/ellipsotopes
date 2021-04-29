function [out,value] = isempty(E,flag_compute_value)
% [out,value] = isempty(E)
%
% This function checks if the ellipsotope E is empty by solving a convex
% program. The output is TRUE if the ellipsotope is EMPTY. The "value"
% output is the cost of the convex program, which we set to 0 in the case
% of an unconstrained ellipsotope.
%
% See also: ellipsotope.contains
%
% Authors: Shreyas Kousik
% Created: 27 Apr 2021
% Updated: 28 Apr 2021

    % get etoproperties
    [p_norm,~,~,A,b,I] = E.get_properties() ;

    if isempty(A)
        % no constraints = not empty
        out = false ;
        value = 0 ;
    else       
        % create initial guess
        x_0 = pinv(A)*b ;
        
        % check if x_0 is actually feasible to the constraints
        if vecnorm(A*x_0 - b) < 1e-10
            out = true ;
            value = inf ;
            warning('The ellipsotope has degenerate constraints!')
        else
            % set up program cost
            cost = @(x) E.cost_for_emptiness_check(x,p_norm,I) ;
            
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
end