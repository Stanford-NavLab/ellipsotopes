function [idxs,log_idxs] = identify_component_ellipsotopes(E)
% idxs = identify_component_ellipsotopes(E)
%
% Find all component ellipsotopes a given ellipsotope. The output indices
% correspond to the etope's index set.
%
% Authors: Shreyas Kousik
% Created: 13 July 2021
% Updated: nope

    % get properties
    [~,~,~,A,~,I] = get_properties(E) ;
    n_I = length(I) ;

    % initialize second output
    log_idxs = true(1,n_I) ;
    
    % if constraint set is empty then every index subset is a component etope
    if isempty(A)
        idxs = 1:n_I ;
    else
        % initialize output
        idxs = [] ;

        % iterate through index set and check for component etopes
        for idx = 1:n_I
            % get current indices
            J = I{idx} ;

            % get all constraints with nonzero values at the current indices
            A_row_log = any(A(:,J) ~= 0,2) ;

            % check if all other entires in the current constraints are zeros
            A_temp = A ;
            A_temp(:,J) = [] ;
            if all(all(A_temp(A_row_log,:) == 0,2))
                idxs = [idxs, idx] ;
            else
                log_idxs(idx) = false ;
            end
        end
    end
end