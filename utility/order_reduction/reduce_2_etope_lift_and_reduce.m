function E_rdc = reduce_2_etope_lift_and_reduce(E,n_rdc)
% E_rdc = reduce_2_etope_lift_and_reduce(E,n_rdc)
%
% Lift a 2-ellipsotope then reduce based on the lifted component
% ellipsotopes. The quantity n_rdc is the number of generators to remove.
%
% See also: reduce_2_etope.m, reduce_2_etope_to_minimal_exact_rep.m
%
% Authors: Shreyas Kousik
% Created: 14 Mar 2022
% Updated: --

    %% setup
    % sanity check
    if E.p_norm ~= 2
        error('This function only works for 2-ellipsotopes!')
    end

    % get... properties... zoinks!
    [~,c,G,A,b,I,n_dim,n_gen,~,n_I] = get_properties(E) ;
    
    % get desired number of output generators
    n_des = n_gen - n_rdc ;

    %% lift
    % lift
    Gl = [G ; A] ;
    cl = [c ; -b] ;

    % save volumes
    vols = nan(1,n_I) ;

    %% reduce to minimal rep
    % NOTE this is different from the exact lift-and-reduce in
    % reduce_2_etope_to_minimal_exact_rep.m
    G_rdc = {} ;
    
    n_gen_remaining = 0 ;
    for idx = 1:n_I
        J = I{idx} ;
        G_idx = Gl(:,J) ;
        [nr,nc] = size(G_idx) ;
        if nc > nr
            G_idx = reduce_2_etope_generator_matrix(G_idx) ;

            % approximate volume of each reduced 'tope
            vols(idx) = ellipsoid_volume_from_generator_matrix(G_idx) ;
        end
        G_rdc = [G_rdc, {G_idx}] ;
        n_gen_remaining = n_gen_remaining + size(G_idx,2) ;
    end

    %% reduce by using MVOEs if possible
    while sum(~isnan(vols)) > 2
        % sort volumes from smallest to largest
        [~,sort_idxs] = sort(vols,'ascend') ;

        %% NOTE TO SELF
        % need to put code here to handle case when many have the same volume
        %% NOTE TO SELF

        % set up to save the reduced 'topes
        n_G = length(G_rdc) ;
        if n_G > 2
            % get the two smallest topes
            G_i = G_rdc{sort_idxs(1)} ;
            G_j = G_rdc{sort_idxs(2)} ;

            % combine them into one
            G_ij = make_MVOE_generator_matrix(G_i,G_j) ;

            % update G_rdc
            sort_idxs([1 2]) = [] ; % delete 'em
            G_rdc = [{G_ij}, G_rdc(sort_idxs)] ;
            vols = [ellipsoid_volume_from_generator_matrix(G_ij), vols(sort_idxs)] ;
        end
    end
    n_gen_remaining = size(cell2mat(G_rdc),2) ;

    %% create output ellipsotope
    % create new index set
    I_rdc = {} ;
    i_last = 0 ;
    for idx = 1:length(G_rdc)
        n_idx = size(G_rdc{idx},2) ;
        I_rdc = [I_rdc, {(1:n_idx) + i_last}] ;
        i_last = get_max_index(I_rdc) ;
    end

    % return reduced generator matrix to array format
    G_rdc = cell2mat(G_rdc) ;

    % make new 'tope
    c_rdc = cl(1:n_dim) ;
    A_rdc = G_rdc((n_dim+1):end,:) ;
    b_rdc = -cl((n_dim+1):end) ;
    G_rdc = G_rdc(1:n_dim,:) ;
    E_rdc = ellipsotope(2,c_rdc,G_rdc,A_rdc,b_rdc,I_rdc) ;
end