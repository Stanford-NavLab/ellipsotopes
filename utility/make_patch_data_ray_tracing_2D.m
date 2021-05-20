function [F,V] = make_patch_data_ray_tracing_2D(p,c,G,A,b,I,n_P)
% [F,V] = make_patch_data_ray_tracing_2D(p,c,G,A,b,I,n_P)
%
% This creates faces and vertices for plotting a 2-D ellipsotope by tracing
% rays outwards from the ellipsotope center to its boundary. The optional
% input argument n_P is the number of rays to trace, by default 200.
%
% Authors: Shreyas Kousik and Adam Dai
% Created: 20 May 2021

    % set default n_P
    if nargin < 7
        n_P = 200 ;
    end

    % get sizes of things
    [n_dim,n_gen] = size(G) ;
    n_con = size(A,1) ;

    % sanity check
    if n_dim ~= 2
        error('This function only works for 2-D ellipsotopes!')
    end

    % initial ray direction
    g = 2*rand(2,1) - 1 ;

    % create small rotation
    th = 2*pi / n_P ;
    R = rotation_matrix_2D(th) ;

    % save lambdas and g vectors
    g_all = nan(2,n_P-1) ;
    lm_all = nan(1,n_P) ;
    
    % initial guess
    if ~isempty(A)
        x_opt = [1 ; pinv(A)*b ] ;
    else
        x_opt = [1 ; zeros(n_gen,1)] ;
    end
    
    % set up optimization options
    options = optimoptions('fmincon','Display','off',...
        'SpecifyObjectiveGradient',true,...
        'SpecifyConstraintGradient',true) ;

    % perform ray tracing
    for g_idx = 1:n_P
        % rotate g
        g = R*g ;
        g_all(:,g_idx) = g ;

        % set up program
        x_0 = x_opt ; % initial guess from previous solution
        cost = @(x) ray_cost(x) ;
        cons = @(x) ray_nonlcon(x,p,I) ;

        A_eq = [zeros(n_con,1), A ;
                -g, G] ;
        b_eq = [b ; zeros(2,1)] ;

        % call fmincon
        x_opt = fmincon(cost,x_0,[],[],A_eq,b_eq,[],[],cons,options) ;
        lm_all(g_idx) = x_opt(1) ;
    end

    % construct the boundary
    B = repmat(c,1,n_P) + g_all.*repmat(lm_all,2,1) ;

    % create output
    F = [1:n_P,1] ;
    V = B' ;
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