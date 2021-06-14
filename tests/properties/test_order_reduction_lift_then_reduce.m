%% description
clear ; clc ;
%% user parameters
% rng
rng(0)

% etope specs
p_norm = 2 ;
n_dim = 2 ;
n_gen = 10 ;
n_con = 0 ;
n_I = 3 ;

% number of generators to reduce; note this must be > n_dim
n_rdc = 3 ;

%% automated from here
% make random etope
[E,c,G,A,b,I] = make_random_ellipsotope(p_norm,n_dim,n_gen,n_con,n_I) ;

% get n_rdc smallest generators
d_G = vecnorm(G,2) ;
[~,idxs_smallest] = sort(d_G,'ascend') ;

% pop the smallest generators
I_new = [] ;

for idx_J = 1:n_I
    J = I{idx_J} ;
    flag_idx_found = false ;
    
    % check if any of the smallest generators' indices are in J
    for idx_rdc = idxs_smallest(1:n_rdc)
        J_log = (J == idx_rdc) ;
        
        if any(J_log)
            flag_idx_found = true ;
            break
        end
    end
    
    % if the index was found, pop it
    if flag_idx_found
        I_temp = {J(J_log), J(~J_log)} ;
    else
        I_temp = {J} ;
    end    
    
    % assemble the new index set
    I_new = [I_new, I_temp] ;
end

n_I_new = length(I_new) ;

% create new etope with popped gens
E_new = ellipsotope(p_norm,c,G,A,b,I_new) ;

%% reduce popped generators
% get all generators that are "solo"
solo_log = cellfun(@(J) length(J) == 1,I_new) ;
idxs_solo = 1:n_I_new ;
idxs_solo = idxs_solo(solo_log) ;

error('left off here!')

%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ;

plot(E)
plot(E_new,'facecolor','g','edgecolor','g')