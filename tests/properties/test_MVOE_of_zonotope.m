%% description
% This script implements the outer-approximation of a zonotope with an
% ellipsoid using the method in [1].
%
% NOTE! This script requires YALMIP and CORA.
%
% [1] Gaßmann, V. and Althoff, M., 2020, July. Scalable Zonotope-Ellipsoid
%     Conversions using the Euclidean Zonotope Norm. In 2020 American
%     Control Conference (ACC) (pp. 4715-4721). IEEE.
clear ; clc

%% user parameters
% rng seed
rng(0) ;

% zonotope number of generators
n_dim = 2 ;
n_gen = 5 ;

%% automated from here
% make a random 2-D zonotope
G = 2*rand(n_dim,n_gen) - 1 ;

% % example generator matrix from [1,(21)]
% G = [1 -2 2 0 3 1 0 ;
%     0 0 -1 -2 -2 -1 0 ;
%     -2 -1 0 0 -2 1 0 ;
%     1 -1 -1 1 -4 0 5 ;
%     -2 1 0 0 1 0 -3] ;

% [n_dim,n_gen] = size(G) ;
c = zeros(n_dim,1) ;

% compute E_0^(-1/2)
E_0 = n_gen*(G*G') ;
E_0_inv_sqrt = E_0^(-1/2) ;
G_t = E_0_inv_sqrt*G ;

% solve [1,(19)] with YALMIP
lm = sdpvar(n_gen,1) ; % decision variable
options = sdpsettings('verbose',0) ; 
obj = sum(lm) ;
cons = [diag(lm) - G_t'*G_t >= 0, lm >= 0] ;
sol = optimize(cons,obj) ;
r_hat = sum(value(lm)) ;

% compute enclosing ellipsoid
E_enc = r_hat * E_0 ;

%% plotting setup
Z = zonotope(c,G) ;

E = ellipsoid(E_enc,c) ;


%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ;
plot(Z) ;
plot(E) ;