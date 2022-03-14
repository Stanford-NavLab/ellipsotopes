function [out,value,x_feas] = isempty(E,flag_compute_value,method)
% [out,value,x_feas] = E.isempty()
%
% This function checks if the ellipsotope E is empty by solving a convex
% program. The output is TRUE if the ellipsotope is EMPTY. The "value"
% output is the cost of the convex program, which we set to 0 in the case
% of an unconstrained ellipsotope. The x_feas output is a feasible point in
% the ellipsotope.
%
% See also: ellipsotope.contains
%
% Authors: Shreyas Kousik
% Created: 27 Apr 2021
% Updated: 15 Feb 2022 (tried out sqp instead of interior point)

    if nargin < 2
        flag_compute_value = false ;
    end

    if nargin < 3
        method = 'feasibility' ;
    end
    
    if nargout == 3
        flag_compute_value = true ;
    end

    % get etoproperties
    [p_norm,c,G,A,b,I] = E.get_properties() ;

    if isempty(A)
        % no constraints = not empty
        out = false ;
        value = [] ;
        x_feas = c ;
    else
        % create initial guess
        x_0 = pinv(A)*b ;

        % check if x_0 is actually feasible to the constraints
        if vecnorm(A*x_0 - b) > 1e-7
            out = true ;
            value = [] ;
            x_feas = [] ;
            warning('The ellipsotope has degenerate constraints!')
        else
            % test if initial guess cost is < 1 so we can bail out early
            value = E.cost_for_emptiness_check(x_0,p_norm,I,false) ;
            bail_out = value <= 1 ;

            if flag_compute_value || (~bail_out)
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
                            'Algorithm','sqp',...
                            'CheckGradients',false,... % useful to set to true sometimes...
                            'SpecifyObjectiveGradient',true) ;

                        % run optimization, woo!
                        try
                            [x_opt,value] = fmincon(cost,x_0,[],[],A,b,[],[],[],options) ;
                            % output boolean
                            out = value > 1 ;
                            x_feas = c + G*x_opt ;
                        catch
                            out = false ;
                            x_feas = [] ;
                        end
%                         end

                    case 'feasibility'
                        % create cost and nonlcon
                        cost = @(x) E.cost_for_emptiness_check_feas(x,A,b) ;
                        cons = @(x) E.nonlcon_for_emptiness_check_feas(x,p_norm,I) ;

                        % fmincon options
                        options = optimoptions('fmincon','Display','off',...
                            'Algorithm','sqp',...
                            'SpecifyObjectiveGradient',true,...
                            'SpecifyConstraintGradient',true) ;

                        % run optimization
                        [~,value] = fmincon(cost,x_0,[],[],[],[],[],[],cons,options) ;

                        out = abs(value) > 1e-10 ;
                        x_feas = [] ;
                end
            else
                out = false ;
                x_feas = x_0 ;
            end
        end
    end
end