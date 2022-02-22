%% user parameters
% tope specs
n_gen = 5 ;
n_con = 1 ;

% number of rays for plotting
n_ray = 100 ;

%% automated from here
% make a random ellipsotope
E = make_random_ellipsotope(2,2,n_gen,n_con) ;

% get properties
[p,c,G,A,b,I,n_dim,n_gen,n_con,n_I] = E.get_properties() ;

% create an initial point inside the 'tope
p_0 = E.analytic_center ;

% create ray directions
D = rand(n_dim,n_ray) ;

%

%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ;
plot(E)

%% helper functions
function [c,gc] = ray_cost(x,n_gen,n_ray)
    c = -sum(x(1:n_ray)) ;
    gc = [-ones(1,n_ray), zeros(1,n_ray*n_gen)] ;
end

function [c,ceq,gc,gceq] = ray_nonlcon(x,p_norm,I,n_gen,n_ray)
error('left off here!')
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