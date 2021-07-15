function E_out = minkowski_sum_etope_cell(E_cell_in)
% E_out = minkowski_sum_etope_cell(E_cell_in)
%
% Does what it says on the box.
%
% Authors: Shreyas Kousik
% Created: 14 July 2021
% Updated: -

E_out = E_cell_in{1} ;

if length(E_cell_in) > 1
    for idx = 2:length(E_cell_in)
        E_out = E_out + E_cell_in{idx} ;
    end
end