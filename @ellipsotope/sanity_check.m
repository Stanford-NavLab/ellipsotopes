function sanity_check(E)
% sanity_check(E)
%
% Cleans up dimensions/sizes of things and throws an error if there is a
% mismatch.
%
% Authors: Shreyas Kousik
% Created: 16 Mar 2022
% Updated: -

    [~,c,G,A,b,I] = get_properties(E) ;

    n_dim_c = length(c) ;
    [n_dim_G,n_gen_G] = size(G) ;
    n_con_b = length(b) ;
    [n_con_A,n_gen_A] = size(A) ;

    if n_dim_c ~= n_dim_G
        error('Center and generator have different dimensions!')
    end
    
    if n_con_b ~= n_con_A
        error('Constraint matrices are mismatched!')
    end
    
    if (~isempty(A)) && (n_gen_G ~= n_gen_A)
        error('Generator and constraint matrices are mismatched!')
    end
    
    check_index_set_validity(I,G) ;
    
    % set properties
    E.dimension = n_dim_c ;
    E.n_cons = n_con_b ;
end