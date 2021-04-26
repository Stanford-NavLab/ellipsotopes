function [out,value] = contains(E,p_test)
% [out,value] = E.contains(p)
%
% Check if the point p is inside the ellipsotope E, meaning the output
% value is <= 1.
%
% Authors: Shreyas Kousik
% Created: 26 Apr 2021
% Updated: nah

    % get etoproperties
    [p_norm,c,G,A,b,I] = E.get_properties() ;
    
    % compute point containtment constraint matrices
    A_eq = [G ; A] ;
    b_eq = [p_test - c ; b] ;

    % set up program cost
    cost = @(x) point_cointainment_cost(x,p_norm,I) ;
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

%% helper functions
function [c,gc] = point_cointainment_cost(x,p_norm,I)
    n_I = length(I) ;
    x_vals = nan(n_I,1) ;

    for idx = 1:n_I
        x_vals(idx) = sum(x(I{idx}).^p_norm) ;
    end

    % compute cost
    [c,c_idx] = max(x_vals) ;

    % compute cost (sub)gradient
    gc = zeros(1,length(x)) ;
    gc(I{c_idx}) = p_norm*(x(I{c_idx}).^(p_norm-1)) ;
end