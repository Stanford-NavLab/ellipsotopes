function I = combine_indices(I_1,I_2)
% I = combine_indices(I_1,I_2)
% 
% Combine two index sets into one.
%
% Authors: Adam Dai and Shreyas Kousik
% Created: who knows!
% Updated: 14 July 2021

    m1 = length(I_1) ;
    m2 = length(I_2) ; 
    I_1_end = I_1{end}(end) ; 
    I = I_1 ;
    for idx = 1:m2
       I{m1+idx} = I_2{idx} + I_1_end ; 
    end
end

