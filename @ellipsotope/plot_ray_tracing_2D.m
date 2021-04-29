function h_E = plot_ray_tracing_2D(E,n_g_test,varargin)
% h = plot_ray_tracing_2D(n,patch_args_in)
%
% Compute the boundary of the ellipsotope E via ray tracing for plotting,
% using n rays (default is 100). This function creates a patch plot handle,
% so you can pass in all the arguments you would pass to patch. Woo!
%
% Authors: Shreyas Kousik
% Created: 28 Apr 2021
% Updated: nup

    %% setup
    % get etoproperties
    [p_norm,c,G,A,b,I] = E.get_properties() ;
    
    % get sizes of things
    n_dim = size(G) ;
    n_con = size(A,1) ;

    % sanity check
    if n_dim ~= 2
        error('This function only works for 2-D ellipsotopes!')
    end

    % initial ray direction
    g = 2*rand(2,1) - 1 ;

    % create small rotation
    th = 2*pi / n_g_test ;
    R = rotation_matrix_2D(th) ;

    % save lambdas and g vectors
    g_all = nan(2,n_g_test-1) ;
    lm_all = nan(1,n_g_test) ;
    
    % initial guess
    x_opt = [1 ; pinv(A)*b ] ;
    
    % set up optimization options
    options = optimoptions('fmincon','Display','off',...
        'SpecifyObjectiveGradient',true,...
        'SpecifyConstraintGradient',true) ;

    %% ray tracing
    for g_idx = 1:n_g_test
        % rotate g
        g = R*g ;
        g_all(:,g_idx) = g ;

        % set up program
        x_0 = x_opt ; % initial guess from previous solution
        cost = @(x) ray_cost(x) ;
        cons = @(x) ray_nonlcon(x,p_norm,I) ;
        A_eq = [zeros(n_con,1), A ;
                -g, G] ;
        b_eq = [b ; zeros(2,1)] ;

        % call fmincon
        x_opt = fmincon(cost,x_0,[],[],A_eq,b_eq,[],[],cons,options) ;
        lm_all(g_idx) = x_opt(1) ;
    end

    % construct the boundary
    B = repmat(c,1,n_g_test) + g_all.*repmat(lm_all,2,1) ;

    %% plotting
    % faces and vertices
    F = [1:n_g_test,1] ;
    V = B' ;

    % default options
    plot_options = [{'facecolor','b','linewidth',1.5,'edgecolor','b',...
                    'facealpha',0.1,'edgealpha',0.5}, varargin{:}] ;

    h = patch('faces',F,'vertices',V,plot_options{:}) ;

    if nargout > 1
        h_E = h ;
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