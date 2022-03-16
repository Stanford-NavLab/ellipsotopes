function clean_properties(E)
% clean_properties(E)
%
% Clean up any empty index subsets or constraints, and remove any
% constraints with all zero rows.
%
% Authors: Shreyas Kousik
% Created: 20 July 2021
% Updated: 16 Mar 2022

    %% index set cleanup
    % delete any empty index subsets
    I = E.index_set ;
    I_empty_log = cellfun(@isempty,I) ;
    I = I(~I_empty_log) ;
    E.index_set = I ;

    %% constraint cleanup
    % delete any empty constraints
    A = E.constraint_A ;
    b = E.constraint_b ;
    if isempty(A)
        E.constraint_A = [] ;
        E.constraint_b = [] ;
    else
        A_zero_log = any(A ~= 0,2) ;
        E.constraint_A = A(A_zero_log,:) ;
        E.constraint_b = b(A_zero_log) ;
    end
    
    %% sanity check
    E.sanity_check()
end