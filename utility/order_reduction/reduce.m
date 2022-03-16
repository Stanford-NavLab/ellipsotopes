function E_rdc = reduce(E,n_gen_to_remove,n_con_to_keep,flag_force_reduce)
% E_rdc = reduce(E)
% E_rdc = reduce(E,n_rdc)
%
% Reduce the number of generators in the ellipsotope E. If the second input
% argument is specified, we try to remove that many generators.
%
% Authors: Shreyas Kousik
% Created: in days of yore
% Updated: 16 Mar 2022 (added warning about asking for too much reduction)
    

    % set default number of generators to reduce
    if nargin < 2
        n_gen_to_remove = 1 ;
    end
    
    if nargin < 3
        n_con_to_keep = 3 ; % threeee is a magic number
    end
    
    if nargin < 4
        flag_force_reduce = false ;
    end

    % set default output
    E_rdc = E ;

    % bail out if E is not full-dimensional
    if E.n_generators < E.dimension
        warning('Ellipsotope is not full-dimensional! Not reducing!')
    else
        % make sure n_rdc wouldn't make the ellipsotope smaller than the system
        % dimension
        if (E.n_generators - n_gen_to_remove) < E.n_dimension
            n_gen_to_remove = E.n_generators - E.n_dimension ;
            warning(['Cannot reduce ellipsotope to less than its dimension, ',...
                'meaning ',num2str(E.n_dimension),' generators.'])
        end

        if E.p_norm == 2
            E_rdc = reduce_2_etope(E,n_gen_to_remove,...
                n_con_to_keep,...
                flag_force_reduce) ;

            % cleanup
            E.clean_properties
        else
            warning(['Sorry, we are still implementing order reduction for ',...
                'ellipsotopes with p-norm not equal to 2.'])
    end
end