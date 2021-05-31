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
% Updated: 30 May 2021 (sped up some things)

    if nargin < 2
        flag_compute_value = false ;
    end

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
        if vecnorm(A*x_0 - b) > 1e-10
            out = true ;
            value = inf ;
            warning('The ellipsotope has degenerate constraints!')
        else
            % test if initial feasible guess cost is < 1
            value = E.cost_for_emptiness_check(x_0,p_norm,I,false) ;
            
            if value <= 1 && ~flag_compute_value
                % the initial guess is feasible so the etope is nonempty
                out = false ;
            else
                % run an LP with the constrained zonotope version of the
                % etope to check emptiness conservatively
                [f_cost,A_ineq,b_ineq,A_eq,b_eq] = E.make_con_zono_empty_check_LP(A,b) ;
                z_opt = linprog(f_cost,A_ineq,b_ineq,A_eq,b_eq)  ;
                value = z_opt(end) ;
                
                if value <= 1
                    % the constrained zonotope overapproximation is empty
                    % so the tope is empty
                    out = true ;
                else
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
                end
            end
        end
    end
end