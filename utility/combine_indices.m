function I = combine_indices(I_1,I_2)
% I = combine_indices(I_1,I_2)
% 
% Combine 2 index sets

m1 = length(I_1); m2 = length(I_2); 
I_1_end = I_1{end}(end); 
I = I_1;
for i = 1:m2
   I{m1+i} = I_2{i} + I_1_end; 
end

end

