function I = make_random_index_set(n_gen,n_I)
% I = make_random_index_set(n_gen)
%
% Given a number of generators n_gen, make a random index set I for
% making a random ellipsotope
%
% Authors: Shreyas Kousik
% Created: 28 Apr 2021
% Updated: 14 Mar 2022 (made the output index sets more random)

    % pick a random number of index subsets; we pick the mean and standard
    % deviation to be n_gen/4 to encourage "clumpy" index subsets (that
    % contain more than one number)
    if nargin < 2
        n_I = rand_int(1,n_gen,n_gen/4,n_gen/4) ;
    end
    
    % create full index set
    I = {1:n_gen} ;
    n_I_cur = length(I) ;
    
    while n_I_cur < n_I
        % pick a random subset of the index set
        idx = rand_int(1,n_I_cur) ;
        
        % get the current index subset
        J = I{idx} ;
        
        % pick a random point at which to split the current index subset
        if length(J) > 1
            idx_split = rand_int(1,length(J)-1) ;
            J_1 = J(1:idx_split) ;
            J_2 = J((idx_split+1):end) ;
            
            I(idx) = [] ;
            I = [I, {J_1}, {J_2}] ;
        end
        
        % continue
        n_I_cur = length(I) ;
    end
    
    % reorder the index set correctly, since the previous operations
    % shuffled it around quite a bit, ohoho
    idxs = cellfun(@(x) x(1),I) ;
    [~,sort_idxs] = sort(idxs,'ascend') ;
    I = I(sort_idxs) ;
end

