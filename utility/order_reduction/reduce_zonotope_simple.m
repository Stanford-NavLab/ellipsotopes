function [G,n_gen,flag_success] = reduce_zonotope_simple(G,n_rdc)
% [G,n_gen,flag_success] = reduce_zonotope_simple(G,n_rdc)
%
% Reduce generator matrix G using the method in [1].
%
% Input:       G      - initial generator matrix
%              n_rdc  - desired number of generators to get rid of
%
% Output:      G_rdc  - reduced generator matrix
%              n_gen  - number of generators in the output matrix
%              flag   - status flag (0 if successful)
%
% [1] Combastel, C., 2003, September. A state bounding observer based on
%     zonotopes. In 2003 European Control Conference (ECC) (pp. 2589-2594).
%     IEEE.
%
% Authors: Scott et al.
% Created: shrug
% Updated: 15 Mar 2022 (Shreyas updated to use in ellipsotope code)

    % initialize optimistically
    flag_success = 0 ;

    % check inputs
    if nargin < 2
        n_rdc = 1 ; % might as well try to reduce by one generator
    end

    [n_dim,n_gen] = size(G) ;

    % sanity check
    if (n_rdc <= 0)
        return ;
    end
    
    L = abs(G) ;

    % another sanity check
    if rank(G) == 1
        % there is just one generator lol
        G = sum(L,2) ;
    else
        % number of generators to aggregate into n generators
        n_agg = n_gen - n_rdc + 1 ;
        if (n_agg <= n_dim)
            return
        end

        % order generators by nearness to a scaled unit vector
        norm_diff = (sum(L,1) - max(L,[],1)) ;
        [~, sort_idxs] = sort(norm_diff,'ascend') ;

        % aggregate generators
        G_ordered = G(:,sort_idxs);
        L2 = abs(G_ordered(:,1:n_agg));
        mag = sum(L2,2) ;
        G = [diag(mag) G_ordered(:,n_agg+1:n_gen)]  ;
    end

    % update output number of generators
    [~, n_gen] = size(G);
end