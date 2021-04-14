%% description
% Test plotting several example ellipsotopes using the built-in class function
%
% Authors: Adam Dai
% Created: 14 Apr 2021
% Updated: 
clear ; clc ;
%% basic ellipsotope

p = 2 ;
c = zeros(2,1) ;
G = [2 1; 1 2];

figure(1);
E = ellipsotope(p,c,G);
plot(E);

%% random general ellipsotope

% ellipsotope definition (make is 2-D please)
p = 2 ;
c = zeros(2,1) ;
G = 2*rand(2,4) - 1;
A = [-1 1 -1 1] ;
b = 0.5 ;
I = {[1,2],[3,4]} ;

figure(2);
E = ellipsotope(p,c,G,A,b,I);
plot(E);
