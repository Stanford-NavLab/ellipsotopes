%% description
clear ; clc
%% user parameters
% ellipsotope definitions
c_1 = [0;0] ;
G_1 = [+1 -1 ;
       -1 +0] ;
   
c_2 = [0;0] ;
G_2 = [-1 -1 ;
       -1 +2] ;
   
c_3 = [0;0] ;
G_3 = [-2 -1 ;
       -1 +1] ;

%% automated from here
E_1 = ellipsotope(2,c_1,G_1) ;
E_2 = ellipsotope(2,c_2,G_2) ;
E_3 = ellipsotope(2,c_3,G_3) ;

% minkowski sum
E = E_1 + E_2 + E_3 ;

%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on ;

plot(E_1,'color','r')
plot(E_2,'color','b')
plot(E_3,'color','g')
plot(E,'color',0.5*ones(1,3))