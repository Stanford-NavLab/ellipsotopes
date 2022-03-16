function [E_reorg,idx_Z_I,idx_Z] = identify_component_zonotope(E)
% [E_reorg,idx_Z_I,idx_Z] = identify_component_zonotope(E)
%
% Identify component zonotope (for the purpose of order reduction).
%
% Inputs: E is an ellipsotope
%
% Outputs:
%   E_reorg - E but with the zonotope generators coming last
%
%   idx_Z_I - the first index of the index set corresponding to comp zono
%
%   idx_Z - the first (column) index of the generators of the comp zono
%
% So, E_reorg.index_set(idx_Z_I:end) are the index subsets of the component
% zonotope, and similarly E_reorg.generators(:,idx_Z:end) are the
% generators of the component zonotope.
%
% See also: identify_component_ellipsotopes.m
%
% Authors: Shreyas Kousik
% Created: 14 Mar 2022
% Updated: 15 Mar 2022

    % get index set
    [p,c,G,A,b,I,~,n_gen,~,n_I] = get_properties(E) ;

    % get lengths of index set elements
    n_per_J = cellfun(@(J) length(J),I) ;

    % count how many elements are in each index subset
    log_idxs = n_per_J == 1 ;

    % get the index set elements that correspond to the component zonotope
    idxs_all = 1:n_I ;
    idxs_Z = idxs_all(log_idxs) ;
    idxs_E = idxs_all(~log_idxs) ;

    % get the generator indices that correspond to the component/other tope
    idxs_Z_G = cell2mat(I(idxs_Z)) ;
    idxs_E_G = cell2mat(I(idxs_E)) ;

    % reorder the generator and constraint matrices
    G_reorg = G(:,[idxs_E_G,idxs_Z_G]) ;
    if ~isempty(A)
        A_reorg = A(:,[idxs_E_G,idxs_Z_G]) ;
    else
        A_reorg = [] ;
    end

    % reorder the index set so the component zonotope bits come last
    I_reorg = [] ;
    n_last = 0 ;
    for idx = idxs_E
        J = I{idx} ;
        n_J = length(J) ;
        I_reorg = [I_reorg, {(1:n_J) + n_last}] ;
        n_last = get_max_index(I_reorg) ;
    end
    idx_Z = n_last + 1 ;
    I_reorg = [I_reorg, num2cell(idx_Z:n_gen)] ;
    
    idx_Z_I = find(cellfun(@(J) length(J),I_reorg) == 1,1,'first') ;

    % make new reorganized tope
    E_reorg = ellipsotope(p,c,G_reorg,A_reorg,b,I_reorg) ;
end