function [idxs,log_idxs,E_other,Z_comp] = identify_component_zonotopes(E,flag_reorg)
% [idxs,log_idxs,E_other,Z_comp] = identify_component_zonotopes(E,flag_reorg)
%
% Identify component zonotopes for the purpose of order reduction. The
% output indices correspond to the etope's index set.
%
% See also: identify_component_ellipsotopes.m
%
% Authors: Shreyas Kousik
% Created: 14 Mar 2022
% Updated: 15 Mar 2022

    % get those propzzzz
    [p,c,G,A,b,I] = get_properties(E) ;
    n_I = length(I) ;
    
    if nargin < 2
        flag_reorg = true ;
    end
    
    % get lengths of index set elements
    n_per_J = cellfun(@(J) length(J),I) ;
    
    % count how many elements are
    log_idxs = n_per_J == 1 ;
    
    % get the index set elements that correspond to the component zonotope
    idxs = 1:n_I ;
    idxs = idxs(log_idxs) ;
    
    % extract the component zonotope
    E_other = {} ;
    Z_comp = {} ;
    
    error('Still gotta finish this!')
end