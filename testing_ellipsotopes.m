%% user parameters
% generator matrix (this should always have three columns)
G = [1 0 2;
     2 1 1] ;
G = pinv([0.6056   -0.2487 ;
   -0.2487    0.5434]) ;
%  % constraints
A = [0 -0.1 0.1] ;
b = 0 ;
tol = 1e-2 ;

% norm to consider
norm_p = 2 ;

% number of points to use for generating the plot
n_P = 30 ;
 
%% automated from here
% get number of generators
n_G = size(G,2) ;

B = make_grid(repmat([-1,1],1,n_G),n_P*ones(1,n_G)) ;

% % evaluate which points obey the constraints
C_log = abs(A*B - b) <= tol ;
B = B(:,C_log) ; 
% sum(C_log)

% evaluate which points obey the norm
N_log = vecnorm(B,norm_p) <= 1 ;
B = B(:,N_log) ;

% get all the points
P = G*B ;
K_P = boundary(P') ;

%% comparison to ellipsoid
% taking a guess at the PSD matrix... HUZZAH
G_Q = pinv(G) ;
Q = (G_Q'*G_Q) ;

% get grid of points
X_bounds = [min(P(1,:)),max(P(1,:)),min(P(2,:)),max(P(2,:))] ;
X = make_grid(X_bounds) ;
n_X = size(X,2) ;

% evaluate x * Q * x' for each column of X
V = nan(1,n_X) ;
for idx = 1:n_X
    x = X(:,idx) ;
    V(idx) = x'*Q*x ;
end

V_log = V <= 1 ;

X = X(:,V_log) ;
K_X = boundary(X') ;

%% map of a unit circle
[F,V] = make_circle() ;
V = (G*V')' ;

%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ;

plot_path(P,'r.')
plot_path(X,'g.')
plot_path(P(:,K_P),'b')
plot_path(X(:,K_X),'b')
plot_path(V','k')