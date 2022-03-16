function E_rdc = reduce_constraint(E,n_cons_to_keep)
% E_rdc = reduce_constraint(E,n_cons)
% E_rdc = reduce(E,n_rdc)
%
% Reduce the number of constraints in ellipsotope to desired number
% n_cons_to_keep
%
% Authors: Adam Dai
% Created: 9 Mar 2022
% Updated: 

n_rdc = size(E.constraint_A,1) - n_cons_to_keep;

E_rdc = E;

while n_rdc > 0
    % iteratively eliminate one constraint (and generator)
    E_rdc = reduce_etope_constraint_and_generator(E_rdc);
    n_rdc = n_rdc - 1;

end

