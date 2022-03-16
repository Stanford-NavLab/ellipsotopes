%% description
%
clear ; clc ;
%% user parameters
% rng seed
rng(0)

% etope specs
p = 2 ;
n_dim = 2 ; 
n_gen = 15 ;
n_con = 3 ;
n_I = 8 ; % ensure that at least one J \in I has more than one generator

%% automated from here
% make the tope
[E,c,G,A,b,I] = make_random_ellipsotope(p,n_dim,n_gen,n_con,n_I) ;

%% setup
% lift the tope (just the generators
G_l = [G ; A] ;

% get all poppable generators
L = get_index_set_lengths(I) ;
log_pop = L > 1 ;
idx_pop = cell2mat(I(log_pop)) ;

% compute the lengths of poppable generators of the lifted tope
v = vecnorm(G_l(:,idx_pop)) ;
[~,sort_idxs] = sort(v,'ascend') ;
idx_pop = idx_pop(sort_idxs) ;

%% popping
% set up to store each popped tope
n_pop = length(idx_pop) ;
E_pop = cell(1,n_pop) ;

% iterate through poppable generators and pop each one, one at a time
idx = 1 ;
for idx_gen = idx_pop
    E_pop{idx} = pop_generator(E,idx_gen) ;
    idx = idx + 1 ;
end

%% plotting
figure(1) ; clf ;

subplot(1,n_pop+1,1) ;
axis equal ; grid on ;
plot(E) ;

for idx = 1:n_pop
    subplot(1,n_pop+1,idx+1) ;
    axis equal ; grid on ;
    plot(E_pop{idx}) ;
end
