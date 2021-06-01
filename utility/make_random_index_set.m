function I = make_random_index_set(n_gen,n_I)
% I = make_random_index_set(n_gen)
%
% Given a number of generators n_gen, make a random index set I for
% making a random ellipsotope
%
% Authors: Shreyas Kousik
% Created: 28 Apr 2021
% Updated: 31 May 2021 (added n_I input)

    % pick a random number of index subsets; we pick the mean and standard
    % deviation to be n_gen/4 to encourage "clumpy" index subsets (that
    % contain more than one number)
    if nargin < 2
        n_I = rand_int(1,n_gen,n_gen/4,n_gen/4) ;
    end
    
    % figure out how indices to put in each subset
    n_per_J = floor(n_gen/n_I) ;
    
    % create output
    I = cell(1,n_I) ;
    
    n_last = 1 ;
    
    for idx = 1:n_I
        % create index subset
        n_end = (n_last + n_per_J - 1) ;
        if idx == n_I && n_end < n_gen
            n_end = n_gen ;
        end
        J = n_last:n_end ;
        
        % fill in
        I{idx} = J ;
        
        % update last
        n_last = n_last + n_per_J ;
    end
end