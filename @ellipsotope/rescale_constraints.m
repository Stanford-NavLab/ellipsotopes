function E_orig = rescale_constraints(E,scaling_factor)
% rescale_constraints(E,scaling_factor)
%
% This can help with numerical conditioning, you know!
%
% Authors: Shreyas Kousik
% Created: 14 Mar 2022
% Updated: 16 Mar 2022 (tested out removing tiny tiny terms, doesn't work)

    if nargin < 2
        scaling_factor = 1 ;
    end
    
    A = E.constraint_A ;
    b = E.constraint_b ;
    
    V = scaling_factor.*vecnorm([A,b],2,2) ;
    
    if nargout > 0
        E_orig = copy_ellipsotope(E) ;
    end
    
    A_new = A ./ repmat(V,1,E.n_generators) ;
    b_new = b ./ V ;
    
    E.constraint_A = A_new ;
    E.constraint_b = b_new ;
end