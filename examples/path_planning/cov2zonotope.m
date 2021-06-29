function Z = cov2zonotope(cov,m,n)
% cov2zonotope - Converts a covariance matrix to a zonotope with mSigma
% confidence level
%
% Syntax:  
%    Z = cov2zonotope(cov,m)
%
% Inputs:
%    cov - covariance matrix
%    m - confidence interval
%    n - dimensions
%
% Outputs:
%    Z - zonotope representing mSigma interval of covariance matrix
%
% Example: 
%
% Other m-files required: 
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Akshay Shetty
% Written:      15-June-2019
% Last update:---
% Last revision:---

%------------- BEGIN CODE --------------

eps = sqrt(chi2inv(erf(m/sqrt(2)),n));

dim = size(cov,1);

[V,D] = eig(cov);
D = round(D,4);
G = V*sqrt(D);

newG = eps*G;
Z = zonotope([zeros(dim,1), newG]);

%------------- END OF CODE --------------