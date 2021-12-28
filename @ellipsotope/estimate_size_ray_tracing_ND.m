function [v,P,V,D,p_out,x_out] = estimate_size_ray_tracing_ND(E,D_or_n_ray,p_0,x_0)
% v = estimate_size_ray_tracing_ND(E)
%
% Estimate the size of a (nonempty) ellipsotope by tracing rays outwards
% from a point p_0 until they hit the boundary, then averaging the lengths
% of the rays. So, it's not really a volume per se, but it works.
%
% Authors: Shreyas Kousik
% Created: 23 Dec 2021
% Updated: 28 Dec 2021 (forgot to divide by n_rays lol)

if ~isempty(E)
    [p_norm,c,G,A,b,I,n_dim,n_gen,n_con] = E.get_properties() ;
    
    % get a point inside the 'tope
    if nargin < 3
        [p_0,x_0] = E.analytic_center() ;
    end
    
    % create a bunch of random ray directions
    if nargin < 2
        n_ray = 100 ;
        D = 2.*rand(n_dim,n_ray) - 1 ;
        D = D./vecnorm(D) ;
    else
        if length(D_or_n_ray(:)) == 1
            n_ray = D_or_n_ray ;
            D = 2.*rand(n_dim,n_ray) - 1 ;
            D = D./vecnorm(D) ;
        else
            D = D_or_n_ray ;
            n_ray = size(D,2) ;
        end
    end
    
    % set up to save output data
    P = nan(n_dim,n_ray) ;
    V = nan(1,n_ray) ;
    
    % set up optimization options
    options = optimoptions('fmincon','Display','off',...
        'SpecifyObjectiveGradient',true,...
        'SpecifyConstraintGradient',true) ;
    
    % for each ray direction, trace until the boundary of the tope
    for idx = 1:n_ray
        % set up program
        x_0_idx = [1 ; x_0] ; % initial guess from previous solution
        cost = @(x) ray_cost(x) ;
        cons = @(x) ray_nonlcon(x,p_norm,I) ;
        
        A_eq = [zeros(n_con,1), A ;
            -D(:,idx), G] ;
        b_eq = [b ; p_0 - c] ;
        
        % call fmincon
        x_opt = fmincon(cost,x_0_idx,[],[],A_eq,b_eq,[],[],cons,options) ;
        V(idx) = x_opt(1) ;
        
        % save output point
        P(:,idx) = p_0 + D(:,idx).*V(idx) ;
    end
    
    % output size
    v = sum(V)./n_ray ;
    
    % output starting point and starting coefficient
    p_out = p_0 ;
    x_out = x_0 ;
else
    v = 0 ;
    P = [] ;
    V = [] ;
    D = [] ;
    p_out = [] ;
    x_out = [] ;
end
end

%% helper functions
function [c,gc] = ray_cost(x)
    c = -x(1) ;
    gc = [-1, zeros(1,length(x)-1)] ;
end

function [c,ceq,gc,gceq] = ray_nonlcon(x,p_norm,I)
    % original constraint computation method
    % c = cellfun(@(J) sum(x(J).^p_norm),I) - 1 ;

    % get sizes of things
    n_x = length(x) ;
    n_I = length(I) ;
    
    % preallocate constraint and gradient
    c = zeros(n_I,1) ;
    gc = zeros(n_x,n_I) ;
    
    % evaluate nonlinear constraints on coef (||coef(J)||_p <= 1)
    z = x(2:end) ; % coefs (first dec var is scalar value)
    
    for idx = 1:n_I
        J = I{idx} ;
        c(idx) = sum(z(J).^p_norm) - 1 ;
        gc(J+1,idx) = p_norm*z(J).^(p_norm-1) ;
    end
    
    % equality constraint outputs
    ceq = [] ;
    gceq = [] ;
end