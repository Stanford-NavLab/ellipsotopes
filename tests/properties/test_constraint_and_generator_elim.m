%% description
% This script tests overapproximating an ellipsotope per [1, Prop. 5],
% which allows us to delete a constraint and a generator
%
% [1] Scott, J.K., Raimondo, D.M., Marseglia, G.R. and Braatz, R.D., 2016.
%     Constrained zonotopes: A new tool for set-based estimation and fault
%     detection. Automatica, 69, pp.126-136.
%
% Authors: Shreyas Kousik
% Created: 29 May 2021
clear ; clc
%% user parameters
% rng
rng(0)

% etope specs
p_norm = 2 ;
n_dim = 2 ;
n_gen = 10 ;
n_con = 3 ;

% which generator to eliminate
j_rdc = 3 ;

%% autoamted from here
% make random etope
[E,c,G,A,b,I] = make_random_ellipsotope(p_norm,n_dim,n_gen,n_con) ;

vecnorm(A(1,:)./G)

% solve the first constraint for the coefficient value
E_j1 = zeros(n_gen,n_con) ;
E_j1(j_rdc,1) = 1 ;
a_1j_inv = 1./A(1,j_rdc) ;

% create Gamma and Lambda
Gm = G*E_j1*a_1j_inv ;
Lm = A*E_j1*a_1j_inv ;

% create new etope
c_rdc = c + Gm*b ;
G_rdc = G - Gm*A ;
A_rdc = A - Lm*A ;
b_rdc = b - Lm*b ;

% delete first constraint and j-th generator column
A_rdc = A_rdc(2:end,:) ;
b_rdc = b_rdc(2:end) ;
A_rdc(:,j_rdc) = [] ;

% delete j-th generator
G_rdc(:,j_rdc) = [] ;

% reorganize the index set
n_I = length(I) ;
I_rdc = cell(1,n_I) ;
for idx = 1:n_I
    J = I{idx} ;
    J_log = J == j_rdc ;
    
    if any(J_log)
        J(J_log) = [] ;
    end
    
    J_log_minus = J > j_rdc ;
    J(J_log_minus) = J(J_log_minus) - 1 ;
    
    I_rdc{idx} = J ;
end

E_rdc = ellipsotope(p_norm,c_rdc,G_rdc,A_rdc,b_rdc,I_rdc) ;

%% compute bounds on ball product intersecting affine subspace
% initialize
int_E = [-ones(n_gen,1), ones(n_gen,1)] ;
int_R = [-inf(n_gen,1), inf(n_gen,1)] ;

% iterate Algorithm 1 of [1] NOTE! THIS IS BROKEN RIGHT NOW FOR SOME REASON
for idx_i = 1:n_con
    for idx_j = 1:n_gen
        a_ij = A(idx_i,idx_j) ;
        
        if a_ij ~= 0
            R_j = int_R(idx_j,:) ;
            
            a_ij_inv = 1./a_ij ;
            
            idxs_k = 1:n_gen ;
            idxs_k(idx_j) = [] ;
                        
            R_j_RHS = a_ij_inv*b(idx_i).*[1 1] ;
            for idx_k = idxs_k
                a_ik = A(idx_i,idx_k) ;
                E_k = int_E(idx_k,:) ;
                R_j_RHS = interval_subtract(R_j_RHS, a_ij_inv.*a_ik.*E_k) ;
            end
            
            int_R(idx_j,:) = interval_intersect(R_j, R_j_RHS) ;
            
            int_E(idx_j,:) = interval_intersect(int_E(idx_j,:),R_j) ;
        end
    end
end

% %% test computing approximate Hausdorff distance
% % this applies the method in the appendix of [1]
% 
% % compute Q submatrix
% Q = G'*G + eye(n_gen) ;
% Qinv = inv(Q) ;
% 
% % compute H_hat approx for each j
% for j_rdc = 1:n_gen
%     % create j unit vector
%     e_j = zeros(n_gen + n_gcon,1) ;
%     e_j(j_rdc) = 1 ;
%     
%     % create matrix
% end

%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ;

% plot reduced etope
plot(E_rdc,'facecolor','r','edgecolor','r','linestyle','--','facealpha',0.1)

% plot origetope
plot(E) ;

%% helper functions
function I = interval_intersect(I_1,I_2)
     I_lo = [I_1(:,1), I_2(:,1)] ;
     I_hi = [I_1(:,2), I_2(:,2)] ;
     
     I = [max(I_lo,[],2), min(I_hi,[],2)] ;
end

function I = interval_add(I_1,I_2)
    I = [I_1(:,1) + I_2(:,1), I_1(:,2) + I_2(:,2)] ;
end

function I = interval_subtract(I_1,I_2)
    I = [I_1(:,1) - I_2(:,2), I_1(:,2) - I_2(:,1)] ;
end