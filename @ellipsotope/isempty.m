function [out,value] = isempty(E,flag_compute_value,method)
% [out,value] = E.isempty()
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
% Updated: 1 Jun 2021 (added method from corollary in paper)

    if nargin < 2
        flag_compute_value = false ;
    end

    if nargin < 3
        method = 'feasibility' ;
    end

    % get etoproperties
    [p_norm,~,~,A,b,I] = E.get_properties() ;

    if isempty(A)
        % no constraints = not empty
        out = false ;
        value = [] ;
    else
        % create initial guess
        x_0 = pinv(A)*b ;

        % check if x_0 is actually feasible to the constraints
        if vecnorm(A*x_0 - b) > 1e-10
            out = true ;
            value = inf ;
            warning('The ellipsotope has degenerate constraints!')
        else
            % test if initial guess cost is < 1 (i.e., bail out early)
            value = E.cost_for_emptiness_check(x_0,p_norm,I,false) ;
            out = value > 1 ;

            if flag_compute_value || (~out)
                switch method
                    case 'standard'
%                         % run an LP with the constrained zonotope version of the
%                         % etope to check emptiness conservatively
%                         [out,value] = E.isempty_bounding_zonotope() ;
% 
%                         if ~out
                        % set up program cost
                        cost = @(x) E.cost_for_emptiness_check(x,p_norm,I,~flag_compute_value) ;

                        % set up options
                        options = optimoptions('fmincon','Display','off',...
                            'CheckGradients',false,... % useful to set to true sometimes...
                            'SpecifyObjectiveGradient',true) ;

                        % run optimization, woo!
                        try
                            [~,value] = fmincon(cost,x_0,[],[],A,b,[],[],[],options) ;
                            % output boolean
                            out = value > 1 ;
                        catch
                            out = false ;
                        end
%                         end

                    case 'feasibility'
                        % create cost and nonlcon
                        cost = @(x) E.cost_for_emptiness_check_feas(x,A,b) ;
                        cons = @(x) E.nonlcon_for_emptiness_check_feas(x,p_norm,I) ;

                        % fmincon options
                        options = optimoptions('fmincon','Display','off',...
                            'SpecifyObjectiveGradient',true,...
                            'SpecifyConstraintGradient',true) ;

                        % run optimization
                        [~,value] = fmincon(cost,x_0,[],[],[],[],[],[],cons,options) ;

                        out = abs(value) > 1e-10 ;
                end
            else
                value = [] ;
            end
        end
    end
end