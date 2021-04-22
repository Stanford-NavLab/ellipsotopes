function out = and(E1, E2)
% and(E1, E2)
% E1 & E2
%
% Compute the intersection of 2 ellipsotopes (as a new ellipsotope). 
% Overloads the & operator for 2 ellipsotopes.
% Assumes both ellipsotopes have the same p_norm and are basic.
%
% Authors: Adam Dai 
% Created: 1 Mar 2021 
% Updated: 

% extract properties
[p1,c1,G1,A1,b1,I1] = E1.get_properties;
[p2,c2,G2,A2,b2,I2] = E2.get_properties;
m1 = E1.order; m2 = E2.order;
n1 = size(c1,1); n2 = size(c2,1);

% sanity check
if n1 ~= n2
    error('Ellipsotopes to interesect are of different dimension')
end
n = n1;
if p1 ~= p2
    error(['We do not yet support operations on ellipsotopes with ',...
        'different p-norms!'])
end
p = p1;

c = c1;
G = [G1 zeros(size(G2))];

% both basic
if E1.is_basic() && E2.is_basic()
    A = [G1 -G2];
    b = c2 - c1;
    I = {1:m1,m1+1:m1+m2};
% general case
else
    A = [A1                   zeros(size(A1,1),m2);
         zeros(size(A2,1),m1) A2;
         G1                  -G2];
    b = [b1; b2; c2 - c1];
    I = combine_indices(I1, I2);
end

out = ellipsotope(p,c,G,A,b,I);

end