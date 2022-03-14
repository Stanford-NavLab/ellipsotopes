function I = interval_intersect(I_1,I_2)
     I_lo = [I_1(:,1), I_2(:,1)] ;
     I_hi = [I_1(:,2), I_2(:,2)] ;
     
     I = [max(I_lo,[],2), min(I_hi,[],2)] ;
end