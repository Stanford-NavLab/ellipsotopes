%% description
% This script demonstrates the overapproximation created by "popping" a
% single generator. We test out a bunch of different generators that could
% be popped, and maybe take a guess at which one is the least
% overapproximative when popped...? I don't know, man, this is hard.
%
% Authors: Shreyas Kousik
% Created: 16 Mar 2022
% Updated: no
clear ; clc ;
%% user parameters
% rng seed
rng(0)

% etope specs
p = 2 ;
n_dim = 2 ; 
n_gen = 10 ;
n_con = 2 ;
n_I = 3 ; % ensure that at least one J \in I has more than one generator

% number of generators to pop
n_pop = 4 ;

%% automated from here
% make the tope
[E,c,G,A,b,I] = make_random_ellipsotope(p,n_dim,n_gen,n_con,n_I) ;

%% setup
% lift the tope (just the generators
G_l = [G ; A] ;
% G_l = G ;

% get all poppable generators
L = get_index_set_lengths(I) ;
log_pop = L > 1 ;
idx_pop = cell2mat(I(log_pop)) ;

%% try out different heuristics for which gen to pop
% It seems, for now, like picking the lifted generator of smallest
% 2-norm length is the best one to pop... usually. At least, in 2-D for a
% few different random seeds.

% compute the lengths of poppable generators of the lifted tope
G_pop = G_l(:,idx_pop) ;
v = vecnorm(G_pop) ;
[v,sort_idxs] = sort(v,'ascend') ;

% % order generators by nearness to a scaled unit vector
% v = (sum(G_pop,1) - max(G_pop,[],1)) ;
% [v, sort_idxs] = sort(v,'ascend') ;

% reorganize the popperz
idx_pop = idx_pop(sort_idxs) ;

% number of gens to pop
if n_pop > length(idx_pop)
    warning(['Cannot pop as many generators as was requested! ',...
        'Popping only ',num2str(length(idx_pop)),' instead.'])
    idx_pop = idx_pop(1:n_pop) ;
end

%% popping
% set up to store each popped tope
E_pop = cell(1,n_pop) ;

% iterate through poppable generators and pop each one, one at a time
idx = 1 ;
for idx_gen = idx_pop
    E_pop{idx} = pop_generator(E,idx_gen) ;
    idx = idx + 1 ;
end

%% plotting
figure(1) ; clf ; axis equal ; grid on ; hold on ;

colors = [1 0.5 0 ; 0 1 0.5] ; % create colors

for idx = 1:n_pop
    d = (idx-1)/(n_pop-1) ;
    col = interp1([0;1],colors,d) ;
    plot(E_pop{idx},'color',col,'facealpha',0.01,'linewidth',2) ;
end

plot(E,'linestyle',':') ;

%% plotting cleanup
l = [] ;
for idx = 1:n_pop
    l = [l, {num2str(v(idx))}] ;
end
l = [l, {'orig'}] ;

legend(l{:})

make_plot_pretty() ;