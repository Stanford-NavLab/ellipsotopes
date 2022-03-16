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
    
    % find an invertible subset of generators
    GS = sparse(G);
    [~,~,~,q,~]=lu(GS,'vector');
    
    ParTope = G(:,q(1:n_dim));
    Other_Gen = G;
    Other_Gen(:,q(1:n_dim)) = [];
    
    if cond(ParTope)>1e12
       % there is no invertible set; resort to simple reduction
       warning('No invertible subset of generators! Resorting to simple interval overbound.')
       [G,n_gen,flag_success] = reduce_zonotope_simple(G, n_rdc);
       return ;
    end
    
    %At this stage, n generators are in ParTope, while the rest are in
    %Ordered_Gen. The idea is to remove one generator at a time form
    %Ordered_Gen, add it to ParTope to generate a matrix with n+1
    %generators called ParTopePlus, and apply Chisci's method for optimally
    %enclosing ParTopePlus with a parallelotope. This updates ParTope, and
    %the process is repeated until N_Elim generators have been eliminated.
    
    %The loop below attempts to improve the initial parallelotope ParTope
    %by exchanging some generators with Ordered_Gen before any actual
    %reduction is done. Chisci's theorem is used to check, for each
    %generator g in Ordered_Gen, which generator would be eliminated in the
    %optimal reduction of [ParTope | g]. The generator g for which it is
    %most desirable to keep g and eliminate a generator t of ParTope is
    %taken into ParTope, while t is put back in Ordered_Gen. This is
    %repeated at most n times, are until no swap is desirable.
    for idx=1:n_dim
        [SwapDesirabilityMeasure, CondInv] = linsolve(ParTope,Other_Gen);
        SwapDesirabilityMeasure = abs(SwapDesirabilityMeasure);
        
        %Find maximum element of SwapDesirabilityMeasure
        if (CondInv<10e-12)
            warning('Singular matrix should have already been flagged')
        else 
            [ColMax,row_index] = max(SwapDesirabilityMeasure);
            [Max,col_index] = max(ColMax);
            row_index=row_index(col_index);
        end
        
        if (Max>1)
            %It is desirable to swap a generator of ParTope with one of
            %Ordered_Gen. The generators to swap are given by col_index and
            %row_index
            [Other_Gen(:,col_index),ParTope(:,row_index)]...
                =deal(ParTope(:,row_index),Other_Gen(:,col_index));
            
            %Ordered_Gen(:,end+1)=ParTope(:,row_index);
            %ParTope(:,row_index)=Ordered_Gen(:,col_index);
            %Ordered_Gen(:,col_index)=[];
        else
            %It is not desirable to swap any generators
            break;
        end
        
    end
    
    %The loop below eliminates one generator from Ordered_Gen in each pass
    %by adding it to ParTope and reducing back to a parallelotope.
    for idx = 1:n_rdc
        
        %Which generator in Ordered_Gen should be added to ParTope for
        %reduction? Can argue that t is a good generator if column of
        %ReductionDesirabilityMeasure is either very small, very large, or
        %nearly a scaled unit vector. Simple hueristic below considers only
        %the first two cases.
        [ReductionDesirabilityMeasure, CondInv] = linsolve(ParTope,Other_Gen);
        ReductionDesirabilityMeasure = abs(ReductionDesirabilityMeasure);
        if (CondInv<10e-12)
            warning('Cannot Reduce: Singular matrix')
            flag_success = -1;
            return;
        end
        
        %Find max and min element of ReductionDesirabilityMeasure
        [ColMax,row_index] = max(ReductionDesirabilityMeasure);
        [Max,colmax_index] = max(ColMax);
        [Min,colmin_index] = min(ColMax);
        rowmax_index=row_index(colmax_index);
        rowmin_index=row_index(colmin_index);
        
        % Hueristic for choosing which generator to add
        if ( abs(1/Max)<abs(Min) )
            %Very large
            col_index = colmax_index;
            row_index = rowmax_index;
            Val = Max;
        else
            %very small
            col_index = colmin_index;
            row_index = rowmin_index;
            Val = Min;
        end
        
        %Remove generator and from Ordered_Gen
        tnp1 = Other_Gen(:,col_index);
        Other_Gen(:,col_index)=[];
        
        %Reduce n+1 gnerators in ParTopePlus = [ParTope tnp1] to n using
        %Chisci's parallelotope rule.
        %Which generator of ParTopePlus to eliminate?
        RDM = [ReductionDesirabilityMeasure(:,col_index);1];
        if (Val>1)
            ParTope(:,row_index) = tnp1;
            RDM(n_dim+1) = RDM(row_index);
            RDM(row_index) = 1;
        end
        
        %At this point, the generators in ParTope will be preserved.
        %Need to scale these generators to compensate for the eliminated
        %generator
        %r=1+RDM(1:n)/RDM(n+1);
        ParTope = ParTope*diag(1+RDM(1:n_dim)/RDM(n_dim+1));
        
    end
    
    %Reform zonotope generator matrix
    G = [ParTope Other_Gen ];
    
    %Update ng
    [~, n_gen] = size(G);
end