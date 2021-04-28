%% description
% This script tests maximizing the distance along a ray within an
% ellipsotope to find points on the ellipsotope boundary. This lets us plot
% the boundary of an ellipsotope.
%
% The tradeoff between this method of plotting versus the "normal" method
% seems to be 5 generators.
%
% Authors: Shreyas Kousik
% Created 27 Apr 2021
% Updated: 28 Apr 2021 ;
clear ; clc
%% user parameters
% rng seed
rng(0) ;

% ellipsotope properties
p_norm = 2 ;
n_gen = 6 ;

% number of lines to test
n_g_test = 200 ;

%% automated from here
% make random properties
c = rand(2,1) ;
G = 2*rand(2,n_gen) - 1;
A = rand(1,n_gen) ;
b = rand(1) ;
I = make_random_index_set(n_gen) ;

% initial line direction
g = 2*rand(2,1) - 1 ;

% etope
E = ellipsotope(p_norm,c,G,A,b,I) ;

% start timer for boundary computation
t_bdry = tic ;

% get sizes of things
[n_dim,n_gen] = size(G) ;
n_con = size(A,1) ;

% set up program
x_0 = [1 ; zeros(n_gen,1)] ; % initial guess (lmabda, beta)
cost = @(x) ray_cost(x) ;
cons = @(x) ray_nonlcon(x,p_norm,I) ;
A_eq = [zeros(n_con,1), A ;
        -g, G] ;
b_eq = [b ; zeros(n_dim,1)] ;

% set up optimization options
options = optimoptions('fmincon','Display','off',...
    'SpecifyObjectiveGradient',true,...
    'SpecifyConstraintGradient',true) ;

% call fmincon
[x_opt,f_val] = fmincon(cost,x_0,[],[],A_eq,b_eq,[],[],cons,options) ;
lm = x_opt(1) ;

%% test rotating g
% create small rotation
th = 2*pi / n_g_test ;
R = rotation_matrix_2D(th) ;

% save lambdas and g vectors
g_all = [g, nan(2,n_g_test)] ;
lm_all = [lm, nan(1,n_g_test)] ;

for g_idx = 1:n_g_test
    % rotate g
    g = R*g ;
    g_all(:,g_idx+1) = g ;
    
    % set up program
    x_0 = x_opt ; % initial guess from previous solution
    cost = @(x) ray_cost(x) ;
    cons = @(x) ray_nonlcon(x,p_norm,I) ;
    A_eq = [zeros(n_con,1), A ;
            -g, G] ;
    b_eq = [b ; zeros(n_dim,1)] ;
    
    % call fmincon
    [x_opt,f_val] = fmincon(cost,x_0,[],[],A_eq,b_eq,[],[],cons,options) ;
    lm_all(g_idx+1) = x_opt(1) ;
end
t_bdry = toc(t_bdry) ;

% construct the boundary
B = repmat(c,1,n_g_test+1) + g_all.*lm_all ;

disp(['Time elapsed computing boundary: ',num2str(t_bdry,'%0.2f'),' s'])

%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ;

% plot etope beneath everything else
disp('Plotting ellipsotope using points in coefficient space')
if n_gen <= 6
    t_plot_regular = tic ;
    h_E = plot(E) ;
    t_plot_regular = toc(t_plot_regular) ;
    disp(['Time elapsed with regular plot: ',num2str(t_plot_regular,'%0.2f'),' s'])
else
    h_E = [] ;
end

% plot ray
plot_path([c, c+lm*g],'r-')
h_c = plot(c(1),c(2),'ro','markerfacecolor','r') ;
h_arrow = quiver(c(1),c(2),g(1),g(2),'linewidth',1.5,'color','r','maxheadsize',0.5) ;

% plot boundary
disp('Plotting ellipsotope with ray tracing boundary')
h_bdry = plot_path(B,'r-','linewidth',1.5) ;

% labeling and stuff
xlabel('x_1')
ylabel('x_2')
if ~isempty(h_E)
    legend([h_c, h_E, h_arrow, h_bdry],{'center','ellipsotope','ray','boundary'})
else
    legend([h_c, h_arrow, h_bdry],{'center','ray','boundary'})
end
set_plot_fontsize(15) ;

%% helper function
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