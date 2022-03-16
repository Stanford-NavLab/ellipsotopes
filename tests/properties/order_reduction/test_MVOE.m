%% description
% This script implements a method of computing a minimum-volume outer
% ellipsoid (MVOE) for a Minkowski sum of ellipsoids, per [1].
%
% NOTE! In [1] the ellipsoid is defined by the INVERSE of the shape matrix!
%
% [1] Halder, A., 2018, December. On the parameterized computation of
%     minimum volume outer ellipsoid of Minkowski sum of ellipsoids. In
%     2018 IEEE Conference on Decision and Control (CDC) (pp. 4040-4045).
%     IEEE.
%
% Authors: Shreyas Kousik
% Created: 18 May 2021
% Updated: nah
clear ; clc
%% user parameters
% rng seed
rng(0) ;

% define a pair of ellipsotopes
E_1 = ellipsotope(2,zeros(2,1),2*rand(2)-1) ;
E_2 = ellipsotope(2,zeros(2,1),2*rand(2)-1) ;

% tolerate for iteration
tol = 1e-10 ;

%% automated from here
% get the shape matrices for each ellipsotope
G_1 = E_1.generators ;
Q_1 = inv(pinv(G_1)'*pinv(G_1)) ; % NOTE we invert this matrix!

G_2 = E_2.generators ;
Q_2 = inv(pinv(G_2)'*pinv(G_2)) ; % NOTE we invert this matrix!

%% Halder's method
% get R matrix
R = inv(Q_1)*Q_2 ;

% get lambdas (eigenvalues)
lm = eig(R) ;

% perform iteration until convergence
bt = 0 ;
tol_violation = varphi(bt,lm) ;

while tol_violation > tol
    % iterate eq. (21)
    n = sum(1./(1 + bt.*lm)) ;
    d = sum(lm./(1 + bt.*lm)) ;
    bt = sqrt(n/d) ;
    tol_violation = varphi(bt,lm) ;
end

% construct new ellipsotope shape matrix
Q_MVOE = inv((1 + 1/bt).*Q_1 + (1 + bt).*Q_2) ;
G_MVOE = inv(Q_MVOE^(1/2)) ;

% construct MVOE as an ellipsotope
c_1 = E_1.center ;
c_2 = E_2.center ;

E_MVOE = ellipsotope(2,c_1+c_2,G_MVOE) ;

% get minkowski sum etope
E_sum = E_1 + E_2 ;

%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ;

h_MVOE = plot(E_MVOE,'facecolor','b','edgecolor','b','facealpha',0.05) ;
h_sum = plot(E_sum,'facecolor','k','edgecolor','k','facealpha',0.05) ;
h_1 = plot(E_1,'facecolor','r','edgecolor','r','facealpha',0.05) ;
h_2 = plot(E_2,'facecolor','g','edgecolor','g','facealpha',0.05) ;

legend([h_1,h_2,h_MVOE,h_sum],{'E_1','E_2','MVOE','E_{sum}'})

%% helper functions
function val = varphi(b,lm)
   val = sum((1 - (b.^2).*lm)./(1 + b.*lm)) ;
end