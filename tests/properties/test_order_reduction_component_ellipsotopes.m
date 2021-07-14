%% description
% This script demonstrates reducing an ellipsotope's number of generators
% by identifying component ellipsotopes.
%
% Authors: Shreyas Kousik
% Created: 13 July 2021
% Updated: not yet
clear ; clc ;

%% automated from here
% make ellipsotope
rng(0)
E_1 = make_random_ellipsotope(2,2,4,1,1) ;
E_2 = make_random_ellipsotope(2,2,4,1,1) ;
E_3 = make_random_ellipsotope(2,2,5,3,2) ;
E = E_1 + E_2 + E_3 ;

% get propz
[p,c,G,A,b,I,n_dim,n_gen,n_con,n_I] = get_properties(E) ;

% create copies of propz
G_old = G ;
A_old = A ;
b_old = b ;
I_old = I ;

% identify component ellipsotopes
[idxs,log_idxs] = E.identify_component_ellipsotopes() ;

% iterate through component ellipsotopes and reduce each one
E_comp = cell(1,length(idxs)) ;
for idx = 1:n_I
    if log_idxs(idx)
        % extract the current component etope's info
        J = I{idx} ;
        G_idx = G(:,J) ;
        A_idx = A(:,J) ;
        A_log = all(A_idx ~= 0, 2) ;
        A_idx = A_idx(A_log,:) ;
        b_idx = b(A_log) ;
        n_I_idx = size(G_idx,2) ;
        
        % make and store new etope
        E_idx = ellipsotope(p,zeros(n_dim,1),G_idx,A_idx,b_idx,{1:n_I_idx}) ;
        E_comp{idx} = E_idx ;
        
        % nan out the component etope from E
        G_old(:,J) = nan ;
        A_old(:,J) = nan ;
        A_old(A_log,:) = nan ;
        b_old(A_log) = nan ;
        
        error('LEFT OFF HERE 10:22 PM 13 JULY!')
    end
end
