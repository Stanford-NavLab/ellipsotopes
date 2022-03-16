function E_out = pop_generator(E,idx_gen)
% E = pop_generator(E,idx_gen)
%
% Pop the generator given by the index
%
% Authors: Shreyas Kousik
% Created: 16 Mar 2022
% Updated: no

    % get properties
    [p,c,G,A,b,I,~,~,~,n_I] = get_properties(E) ;
    
    %% picking a generator to pop
    if nargin < 2
        % lift the tope (just the generators)
        G_l = [G ; A] ;

        % get all poppable generators
        L = get_index_set_lengths(I) ;
        log_pop = L > 1 ;
        idx_pop = cell2mat(I(log_pop)) ;

        % compute the lengths of poppable generators of the lifted tope
        v = vecnorm(G_l(:,idx_pop)) ;
        [~,sort_idxs] = sort(v,'ascend') ;
        idx_gen = idx_pop(sort_idxs(1)) ;
    end
    
    %% popping
    % get the index subset containing idx_gen
    idx_J = get_index_set_index_containing_generator(I,idx_gen) ;
    J = I{idx_J} ;
    
    % pop the generator in the index set
    if idx_J < n_I
        I_end = I((idx_J+1):end) ;
    else
        I_end = [] ;
    end
    I = [I(1:idx_J-1), {J(1:end-1)}, {J(end)}, I_end] ;
    
    % swap the popped generator and the last generator in the index set
    G(:,[idx_gen, J(end)]) = G(:,[J(end), idx_gen]) ;
    A(:,[idx_gen, J(end)]) = A(:,[J(end), idx_gen]) ;
    
    % reconstruct the output tope
    E_out = ellipsotope(p,c,G,A,b,I) ;
end