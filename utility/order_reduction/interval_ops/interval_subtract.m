function I = interval_subtract(I_1,I_2)
    I = [I_1(:,1) - I_2(:,2), I_1(:,2) - I_2(:,1)] ;
end
