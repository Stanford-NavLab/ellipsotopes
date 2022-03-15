function E_orig = rescale_constraints(E,scaling_factor)
% rescale_constraints(E,scaling_factor)
%
% This can help with numerical conditioning, you know!
%
% Authors: Shreyas Kousik
% Created: 14 Mar 2022
% Updated: nope

    if nargin < 2
        scaling_factor = 1 ;
    end
    
    A = E.constraint_A ;
    b = E.constraint_b ;
    
    V = scaling_factor.*vecnorm([A,b],2,2) ;
    
    if nargout > 0
        E_orig = copy_ellipsotope(E) ;
    end
    
    E.constraint_A = A ./ repmat(V,1,E.n_generators) ;
    E.constraint_b = b ./ V ;
end