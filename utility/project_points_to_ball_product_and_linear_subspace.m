function [P_out,n_P_out] = project_points_to_ball_product_and_linear_subspace(P_in,...
    p_norm,A,b,I,tol_point_on_bdry)
% P_out = project_points_to_ball_product_and_linear_subspace(P_in,p_norm,A,b,I)
% % P_out = project_points_to_ball_product_and_linear_subspace(...,tol)
%
% This function projects a bunch of input points, P_in, to the boundary of
% the intersection of a ball product (defined by p_norm and I) and a linear
% subspace (defined by Ax = b). If A and b are empty, then we just project
% to the boundary of the ball product.
%
% An optional input is tol, which is the required constraint satisfaction
% tolerance (by default 1e-12); in general, the projected points are around
% a machine epsilon away from the boundary.
%
% Authors: Shreyas Kousik
% Created: 9 Apr 2021
% Updated: 7 Jun 2021

% number of index subsets
n_I = length(I) ;

% dimension and number of points
n_P = size(P_in,2) ;

% set default inputs tolerance for point being on boundary
if nargin < 6
    tol_point_on_bdry = 1e-12 ;
end

% initialize output
P_out = [] ;

if ~isempty(A)
    % project points to linear subspace
    P_in = project_points_to_linear_subspace(P_in,A,b) ;
    
    % get a feasible point on the constraint set (we assume this exists for
    % now)
    t_con = pinv(A)*b ;
    
    for idx_P = 1:n_P
        % get current point
        p_idx = P_in(:,idx_P) ;
        
        % set
        coef_sols = [] ;
        for idx_J = 1:n_I
            % get current index set
            J = I{idx_J} ;
            
            % get nullspace vectors for the current dimensions
            u = p_idx(J,:) ;
            
            % construct polynomial (binomial) coefficients given p-norm, since
            % we're solving sum((coef*u(J) - t(J))^p) = 1 for coef
            pows = 0:p_norm ;
            
            % apply binomial theorem
            coef_u = u.^pows(end:-1:1) ;
            coef_t = t_con(J).^pows ;
            coef_binoms = vec_nchoosek(p_norm) ; % binomial coefficients
            
            % sum the terms to get all the coefficients for the (scalar) coef poly
            coef_poly = coef_binoms.*sum(coef_u.*coef_t,1) ;
            
            % subtract 1 from both sides to get the polynomial in standard form
            coef_poly(:,end) = coef_poly(:,end) - 1 ;
            
            % solve for roots
            coef_sol = roots(coef_poly) ;
            
            % remove imaginary roots
            coef_sol_imag = imag(coef_sol) ;
            coef_sol(coef_sol_imag ~= 0) = [] ;
            
            % stack 'em
            coef_sols = [coef_sols, coef_sol(:)'] ;
        end
        
        % create candidate points
        p_test = p_idx*coef_sols + t_con ;
        
        % test to find which point obeys the norm
        v_idx = vecnorm_ball_product(p_test,p_norm,I) ;
        v_test = abs(v_idx - 1) < tol_point_on_bdry ;
        
        % save the points that are ok
        P_out = [P_out, p_test(:,v_test)] ;
    end
else
    P_out = project_points_to_ball_product(P_in,p_norm,I) ;
end

% this procedure doubles the number of points usually
n_P_out = size(P_out,2) ;