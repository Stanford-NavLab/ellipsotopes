function E_drop = drop(E,n_dim)
% E_drop = drop(E,n_dim)
%
% Drop E back down to n_dim dimensions. Only works if E is an indexed tope.
%
% Authors: Shreyas Kousik
% Created: 15 Mar 2022
% Updated: nah

    if E.is_indexed()
        [p,c,G,A,b,I] = get_properties(E) ;
        
        if ~(isempty(A) && isempty(b))
            warning('This ellipsotope has constraints! Things may get wonky!')
        end
        
        G_drop = G(1:n_dim,:) ;
        c_drop = c(1:n_dim) ;
        A_drop = G((n_dim+1):end,:) ;
        b_drop = -c((n_dim+1):end) ;
        
        E_drop = ellipsotope(p,c_drop,G_drop,A_drop,b_drop,I) ;
    else
        error('You can only drop a lifted ellipsotope! Hey!')
    end
end