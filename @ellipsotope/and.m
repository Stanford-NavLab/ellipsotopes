% intersection (overloads &)
% assumes both inputs have same p_norm and are basic
function out = and(E1, E2)
c = E1.center;
G = [E1.generators zeros(size(E1.generators))];
A = [E1.generators -E2.generators];
b = E2.center - E1.center;
I = {1:E1.order,E1.order+1:E1.order+E2.order};
%             if is_constrained(E1)
%
%             else
%                 A = []; b= [];
%             end
%             if is_general(E)
%             else
%                 I = [];
%             end
out = ellipsotope(E1.p_norm, c, G, A, b, I);
end