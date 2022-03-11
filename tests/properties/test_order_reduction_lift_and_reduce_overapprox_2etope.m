%% description
% This script tests order reduction for a 2-ellipsotope starting from some
% given number of generators and then cutting down to a user-specified
% number.
%
% See also: figure_order_reduction_2_etopes.m,
% test_order_reduction_high_dim.m, reduce_2_etope_to_minimal_exact_rep.m
%
% Authors: Shreyas Kousik
% Created: 10 Mar 2022
% Updated: 10 Mar 2022
clear ; clc ; close all ;
%% user parameters
% rng seed
rng(0)

% original etope properties
n_dim = 2 ; % leave this as 2 for plotting
n_gen = 25 ; % default is 20
n_con = 4 ; % default is 3
n_I = 3 ; % shrug

% reduction parameters
n_rdc = 12 ; % default is 12

%% automated from here
disp('Making random 2-ellipsotope')
% make random ellipsotope
p_norm = 2 ;
[E,c,G,A,b,I] = make_random_ellipsotope(p_norm,n_dim,n_gen,n_con,n_I) ;

%% lift and reduce
disp('Lifting and reducing to get minimal representation')
% lift
Gl = [G ; A] ;
cl = [c ; -b] ;
nl = size(Gl,1) ; % dimension of lifted etope

% save volumes
vols = nan(1,n_I) ;

% reduce by iterating through I
G_rdc = {} ; % NOTE this is different from the exact lift-and-reduce
% I_rdc = {} ;
n_gen_remaining = 0 ;
for idx = 1:n_I
    J = I{idx} ;
    G_idx = Gl(:,J) ;
    [nr,nc] = size(G_idx) ;
    if nc > nr
        G_idx = reduce_2_etope_generator_matrix(G_idx) ;
        % I_idx = (1:nl) + get_max_index(I_rdc) ;
    else
        % I_idx = (1:nc) + get_max_index(I_rdc) ;
    end
    G_rdc = [G_rdc, {G_idx}] ;
    % I_rdc = [I_rdc, {I_idx}] ; % we'll reconstruct this later
    n_gen_remaining = n_gen_remaining + size(G_idx,2) ;
    
    % approximate volume of each reduced 'tope
    vols(idx) = ellipsoid_volume_from_generator_matrix(G_idx) ;
end

%% overapproximate with MVOEs
disp('Computing MVOEs to get reduced rep')
while (n_gen_remaining > n_rdc)
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
    
    n_gen_remaining = size(cell2mat(G_rdc),2) ;
end

%% create tope
disp('Creating new reduced ellipsotope')

% create new index set
I_rdc = {} ;
i_last = 0 ;
for idx = 1:length(G_rdc)
    n_idx = size(G_rdc{idx},2) ;
    I_rdc = [I_rdc, {(1:n_idx) + i_last}] ;
    i_last = n_idx ;
end

% return reduced generator matrix to array format
G_rdc = cell2mat(G_rdc) ;

% make new 'tope
c_rdc = cl(1:n_dim) ;
A_rdc = G_rdc((n_dim+1):end,:) ;
b_rdc = -cl((n_dim+1):end) ;
G_rdc = G_rdc(1:n_dim,:) ;
E_rdc = ellipsotope(2,c_rdc,G_rdc,A_rdc,b_rdc,I_rdc) ;

%% plotting
if n_dim == 2
    disp('Plotting!')
    figure(1) ; clf ; axis equal ; hold on ; grid on ;
    
    plot(E)
    plot(E_rdc,'color','r','linestyle','--','linewidth',3)
    
    make_plot_pretty()
end