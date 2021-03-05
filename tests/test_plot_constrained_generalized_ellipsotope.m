%% user parameters
% center
c = [0; 0];
% generator matrix 
G = [1 0 0 0;
     0 1 0 0] ;
%  % constraints
A = [1 0 -1  0;
     0 1  0 -1] ;
b = [1; 0];
tol = 1e-2 ;
% index set
J = {[1,2],[3,4]};

% norm to consider
norm_p = 2 ;

% number of points to use for generating the plot
n_P = 1000 ;

% sampling bounds
b_samp = [-2,2];
 
%% automated from here
% get number of generators
n_G = size(G,2) ;

% generate hyperplane defined by constraints
N = null(A);
null_dim = size(N,2);

B = make_grid(repmat(b_samp,1,null_dim),n_P*ones(1,null_dim)) ;

B = N*B + linsolve(A,b);

% evaluate which points obey the norm
for i = 1:length(J)
    N_log = vecnorm(B(J{i},:),norm_p) <= 1 ;
    B = B(:,N_log) ;
end

% get all the points
P = c + G*B;
K_P = boundary(P') ;

%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ;

plot_path(P,'r.')
plot_path(P(:,K_P),'b')