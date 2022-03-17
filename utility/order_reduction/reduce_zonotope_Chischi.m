function [G,n_gen,flag_success] = reduce_zonotope_Chischi(G,n_rdc)
% [G,n_gen,flag_success] = reduce_zonotope_Chischi(G,n_rdc)
%
% Input:       G      - initial generator matrix                                  
%              n_rdc  - desired number of generators to reduce
% Output:      G_rdc  - reduced generator matrix
%              n_gen  - number of generators in the output matrix
%              flag   - status flag (0 if successful)           
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

    [n_dim,n_gen] = size(G);

    % sanity check
    if (n_rdc <= 0)
        return ;
    end

    % make sure we're not trying to reduce more than n_dim
    if (n_gen - n_rdc) < n_dim
        n_rdc = n_gen - n_dim ;
        warning(['Cannot reduce as many generators as requested! ',...
            'Reducing only ',num2str(n_rdc),' of them.'])
    end

    %Find an invertible subset of generators
    [Dual,flag,~,~] = gauss_elim_partial_solve([G eye(n_dim)],n_dim,n_gen);
    if flag == -4
        % there is no invertible set; resort to simple reduction
        warning('No invertible subset of generators! Resorting to simple interval overbound.')
        [G,n_gen,flag_success] = reduce_zonotope_simple(G, n_rdc);
        return ;
    else
        if (flag)
            return
        end
    end

    %At this stage, n generators are in ParTope, while the rest are in
    %Ordered_Gen. The idea is to remove one generator at a time form
    %Ordered_Gen, add it to ParTope to generate a matrix with n+1
    %generators called ParTopePlus, and apply Chisci's method for optimally
    %enclosing ParTopePlus with a parallelotope. This updates ParTope, and
    %the process is repeated until N_Elim generators have been eliminated.

    %The loop below eliminates one generator from Ordered_Gen in each pass
    %by adding it to ParTope and reducing back to a parallelotope.
    for i=1:n_rdc

        %Which generator in Ordered_Gen should be added to ParTope for
        %reduction? Can argue that t is a good generator if column of
        %ReductionDesirabilityMeasure is either very small, very large, or
        %nearly a scaled unit vector. Simple hueristic below considers only
        %the first two cases.
        %[ReductionDesirabilityMeasure, CondInv] = linsolve(ParTope,Other_Gen);
        R = abs(Dual(:,n_dim+1:n_gen));
        RDM = prod(ones(size(R))+R) - 1 - sum(R);

        %Find max and min element of ReductionDesirabilityMeasure
        [~,col_index] = min(RDM);

        %Remove generator
        Dual(:,n_dim+col_index)=[];
        n_gen = n_gen-1;

        %Reduce n+1 gnerators in ParTopePlus = [ParTope tnp1] to n using
        %Chisci's parallelotope rule.
        Dual(:,n_dim+1:n_gen+n_dim) = diag(1./(1+R(1:n_dim,col_index)))*Dual(:,n_dim+1:n_gen+n_dim);
    end

    %Reform zonotope generator matrix
    G = Dual(:,n_gen+1:n_gen+n_dim)\Dual(:,1:n_gen);
end


function [A,flag,row_pivots,col_pivots] = gauss_elim_partial_solve(A,NR,NC)
% function:    PartialSolve(A,NR,NC);                                     %
% Description:                                                            %
% Input:       A      - nr-by-nc matrix, nr<nc                            %
%              NR     - index of the last row able to pivot               % 
%              NC     - index of the last column able to pivot            % 
% Output:      A      - nr-by-nc matrix, nr<nc                            %
%              flag   - status flag. Zero if successful.                  %

    %---------------------------------------------------------%
    %Initialize
    FuncID = 'PartialSolve';
    flag=0;
    [nr nc] = size(A);
    col_pivots=1:nc;
    row_pivots=1:nr;
    %---------------------------------------------------------%


    
    %---------------------------------------------------------%
    %Check inputs
    if (nr==0 || nr>nc || NC>nc || NR>nr);
        disp(['Error in ' FuncID]);
        return;
    end;
    %---------------------------------------------------------%
    
    
    
    %---------------------------------------------------------%
    %
    tol = max(nr,nc)*eps('double')*norm(A,'inf');
    for row=1:NR
    
        %Find the best pivot
        [~,xi_elim] = max(abs(A(row:NR,row:NC)),[],2);
        xi_elim=row-1+xi_elim;
        A_Scaled=A;
        A_Norms=ones(NR+1-row,1);
        for i=row:NR
            if ( abs(A(i,xi_elim(i-row+1)))>tol ) 
                A_Scaled(i,:) = A(i,:)/A(i,xi_elim(i-row+1));
                A_Norms(i) = norm(A_Scaled(i,:),1)-1;
            else
                %Everything in row i from 1:NC is zero
                A(i,1:NC)=0;
                A_Scaled(i,1:NC)=0;
                A_Norms(i)=Inf;
            end
        end
        [~,index] = sort(A_Norms(row:NR));
        pivot_row = index(1)+row-1;
        pivot_col = xi_elim(index(1));
        
        %Permute
        A([row pivot_row],:) = A([pivot_row row],:);
        row_pivots([row pivot_row]) = row_pivots([pivot_row row]);
        A(:, [row pivot_col]) = A(:, [pivot_col row]);
        col_pivots([row pivot_col]) = col_pivots([pivot_col row]);
        
        %Eliminate
        if ( abs(A(row,row))>tol)
            A(row,row:nc)=A(row,row:nc)/A(row,row);
            
            %Eliminate downward
            for i=row+1:nr
                if (abs(A(i,row))>tol)
                    A(i,:) = A(i,:)-A(i,row)*A(row,:);
                else
                    A(i,row)=0;
                end
            end 
            
            %Eliminate upward
            for i=row-1:-1:1
                if (abs(A(i,row))>1e-10)
                    A(i,row:nc) = A(i,row:nc)-A(i,row)*A(row,row:nc);
                else
                    A(i,row)=0;
                end
            end
        end
        
    end
    %---------------------------------------------------------%
   
end