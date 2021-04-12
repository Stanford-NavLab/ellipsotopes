function out = is_reduced(E)
d = E.dimension ;
o = E.order ;

out = E.is_basic && (d == o) ;
end