function E = times(E1,E2)
% E = times(E1,E2)
% E = E1 .* E2 
%
% Cartesian product of ellipsotopes
%
% Authors: Shreyas Kousik
% Created: 21 Sept 2021


[p1,c1,G1,A1,b1,I1] = E1.get_properties() ;
[p2,c2,G2,A2,b2,I2] = E2.get_properties() ;

if p1 ~= p2
    warning(['Currently can''t support ellipsotopes with different p norms!',...
        ' Picking the higher p-norm!']) ;
    p = max(p1,p2) ;
else
    p = p1 ;
end

c = [c1 ; c2] ;
G = blkdiag(G1,G2) ;
A = blkdiag(A1,A2) ;
b = [b1 ; b2] ;
I = combine_indices(I1,I2) ;

E = ellipsotope(p,c,G,A,b,I) ;
end