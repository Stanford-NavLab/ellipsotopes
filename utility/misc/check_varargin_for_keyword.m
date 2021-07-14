function [val,args_out] = check_varargin_for_keyword(kw,varargin)
% val = check_varargin_for_keyword(kw,varargin)
%
% Iterate through varargin, extract the corresponding value, and return the
% value "val" and the remaining varargin.
%
% Authors: Shreyas Kousik
% Created: 20 May 2021
% Updated: nah

    n = length(varargin) ;
    val = [] ;
    args_out = varargin ;
    for idx = 1:2:n
        kw_idx = varargin{idx} ;
        if strcmpi(kw,kw_idx)
            val = varargin{idx+1} ;
            args_out(idx:idx+1) = [] ;
            break
        end
    end 
end