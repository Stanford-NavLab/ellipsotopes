function [p,x] = analytic_center(E)
% [p,x] = analytic_center(E)
%
% Get the analytic center of the ellipsotope E using the method from the
% Boyd and Vandenberghe textbook.
%
% Authors: Shreyas Kousik
% Created: 23 Dec 2021
% Updated: nope

% check that it's not empty
if ~isempty(E)
    [p_norm,c,G,A,b,I] = E.get_properties() ;
    
    % set up optimization options
    options = optimoptions('fmincon','Display','off',...
        'SpecifyObjectiveGradient',true) ;
    
    % solve for center with fmincon
    x_0 = pinv(A)*b ; % initial guess
    x = fmincon(@(x) analytic_center_cost_fn(x,p_norm,I),...
        x_0,[],[],A,b,[],[],[],options) ;
    
    % get center in workspace
    p = c + G*x ;
else
    p = [] ;
    x = [] ;
end
end

function [c,gc] = analytic_center_cost_fn(x,p_norm,I)
    n_X = length(x) ;
    n_I = length(I) ;
    h = nan(n_I,1) ;
    gc = nan(1,n_X) ;
    for idx = 1:n_I
        J = I{idx} ;
        h(idx) = sum(x(J).^p_norm) - 1 ;
        gc(J) = -(p_norm*x(J).^(p_norm-1))./h(idx) ;
    end

    c = -sum(log(-h)) ;
end