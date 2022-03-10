function E_rdc = reduce_2_etope_to_minimal_exact_rep(E)
% E_min = reduce_2_etope_to_minimal_exact_rep(E)
%
% Implementing the magical proposition from the paper!
%
% Authors: Shreyas Kousik
% Created: 10 Mar 2022
% Updated: nah

    %% setup and sanity check
    % get properties
    [p,c,G,A,b,I,n_dim,~,~,n_I] = get_properties(E) ;

    if p ~= 2
        error('This function only works for 2-ellipsotopes!')
    end

    %% lift and reduce
    % lift
    Gl = [G ; A] ;
    cl = [c ; -b] ;
    nl = size(Gl,1) ; % dimension of lifted etope

    % reduce by iterating through I
    G_rdc = [] ;
    I_rdc = {} ;
    for idx = 1:n_I
        J = I{idx} ;
        G_idx = Gl(:,J) ;
        [nr,nc] = size(G_idx) ;
        if nc > nr
            G_idx = reduce_2_etope_generator_matrix(G_idx) ;
            I_idx = (1:nl) + get_max_index(I_rdc) ;
        else
            I_idx = (1:nc) + get_max_index(I_rdc) ;
        end
        G_rdc = [G_rdc, G_idx] ; 
        I_rdc = [I_rdc, {I_idx}] ;
    end

    % make new 'tope
    c_rdc = cl(1:n_dim) ;
    A_rdc = G_rdc((n_dim+1):end,:) ;
    b_rdc = -cl((n_dim+1):end) ;
    G_rdc = G_rdc(1:n_dim,:) ;
    E_rdc = ellipsotope(2,c_rdc,G_rdc,A_rdc,b_rdc,I_rdc) ;
end