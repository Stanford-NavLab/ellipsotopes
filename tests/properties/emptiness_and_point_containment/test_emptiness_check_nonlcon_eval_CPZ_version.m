%% description
% This script tests the CPZ representation of the ellipsotope ball product
% constraints, and compares the solve speed to the current for-loop version
% in @ellipsotope.nonlcon_for_emptiness_check_feas
%
% Authors: Shreyas Kousik
% Created: 21 June 2021
% Updated: 23 June 2021
clear ; clc
%% user parameters
% rng seed
% rng(102)

% ellipsotope specs
p_norm = 2 ;
n_dim = 3 ;
n_gen = 30 ;
n_con = 10 ;
n_I = 20 ;

%% automated from here
% make the tope
[E,c,G,A,b,I] = make_random_ellipsotope(p_norm,n_dim,n_gen,n_con,n_I) ;

% create powers and coefficient matrices
P = zeros(n_I,n_gen) ;
GP = zeros(n_I,n_gen) ;
C = zeros(n_I,n_gen) ;
for idx = 1:n_I
    J = I{idx} ;
    P(idx,J) = p_norm ;
    GP(idx,J) = p_norm - 1 ;
    C(idx,J) = 1 ;
end

%% evaluate constraint both ways
x = 2*rand(n_gen,1) ;

[c_1,~,gc_1,~] = nonlcon_for_emptiness_check_feas(E,x,p_norm,I) ;
[c_2,~,gc_2,~] = nonlcon_CPZ_version(x,C,P,GP) ;

% test validity
if all(c_1 == c_2)
    disp('Constraints match!')
end

if all(gc_1(:) == gc_2(:))
    disp('Constraint gradients match!')
end

%% time the two versions to figure out the "crossover point"
disp('Timing the two nonlcon eval versions')

n_gens = 1:30:800 ;
n_topes_per_gen = 10 ;
n_n_gens = length(n_gens) ;
t_data = nan(n_n_gens,2) ;
t_idx = 1 ;

for n_gen = n_gens
    t_1 = nan(1,n_topes_per_gen) ;
    t_2 = nan(1,n_topes_per_gen) ;
    for idx_tope = 1:n_topes_per_gen
        % make the tope
        [E,c,G,A,b,I] = make_random_ellipsotope(p_norm,n_dim,n_gen,n_con,n_I) ;

        % create powers and coefficient matrices
        P = zeros(n_I,n_gen) ;
        GP = zeros(n_I,n_gen) ;
        C = zeros(n_I,n_gen) ;
        for idx = 1:n_I
            J = I{idx} ;
            P(idx,J) = p_norm ;
            GP(idx,J) = p_norm - 1 ;
            C(idx,J) = 1 ;
        end
        x = 2*rand(n_gen,1) ;
        t_1(idx_tope) = timeit(@() nonlcon_for_emptiness_check_feas(E,x,p_norm,I)) ;
        t_2(idx_tope) = timeit(@() nonlcon_CPZ_version(x,C,P,GP)) ;
    end
    t_data(t_idx ,1) = mean(t_1) ;
    t_data(t_idx ,2) = mean(t_2) ;
    t_idx  = t_idx  + 1 ;
    disp([num2str(100*(t_idx-1)/n_n_gens,'%0.2f'),'% complete'])
end

%%
fh = figure(1) ; clf ; hold on ;
plot(n_gens,t_data)
legend('method 1','method 2','location','best')
make_plot_pretty()

%% helper functions
function [c,ceq,gc,gceq] = nonlcon_CPZ_version(x,C,P,GP)

    % setup
    n_con = size(P,1) ;
    X = repmat(x(:)',n_con,1) ;
    
    % computee constraint and gradient
    c = sum(C.*(X.^P),2) - 1;
    gc = (C.*P.*(X.^GP))' ;
    
    % equality constraints
    ceq = [] ;
    gceq = [] ;
end