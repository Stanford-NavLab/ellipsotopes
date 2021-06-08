function out = plus(s1, s2)
% plus(E1, E2)
% E1 + E2
%
% Computes the Minkowski sum of 2 ellipsotopes (as a new ellipsotope). 
% Overloads the + operator for 2 ellipstopes.
%
% Authors: Adam Dai 
% Created: 1 Mar 2021 
% Updated: 25 Mar 2021

% determine ellipsotope
if isa(s1,'ellipsotope')
    % ellipsotope + vector
    if isnumeric(s2)
        c = s1.center + s2;
        out = ellipsotope(s1.p_norm,c,s1.generators,s1.constraint_A,s1.constraint_b,s1.index_set);
        return
    % ellipsotope + ellipsotope
    elseif isa(s2,'ellipsotope')
        % both basic
        if is_basic(s1) && is_basic(s2)
            c = s1.center + s2.center;
            G = [s1.generators s2.generators];
            I = {1:s1.order,s1.order+1:s1.order+s2.order};
            out = ellipsotope(s1.p_norm,c,G,[],[],I);
            return
        % general case (constrained and indexed)
        else
            c = s1.center + s2.center;
            G = [s1.generators s2.generators];
            I = combine_indices(s1.index_set, s2.index_set);
            % handle constrained/unconstrained cases
            if s1.is_constrained()
                A1_con = s1.constraint_A;
            else
                A1_con = zeros(0,s1.order);
            end
            if s2.is_constrained()
                A2_con = s2.constraint_A;
            else
                A2_con = zeros(0,s2.order);
            end
            A = blkdiag(A1_con, A2_con);
            b = [s1.constraint_b; s2.constraint_b];
            out = ellipsotope(s1.p_norm,c,G,A,b,I);
            return
        end
    end
% vector + ellipsotope
elseif isnumeric(s1)
    c = s2.center + s1;
    out = ellipsotope(s2.p_norm,c,s2.generators,s2.constraint_A,s2.constraint_b,s2.index_set);
    return
end

end