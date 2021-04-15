%% description
% Test plotting several example ellipsotopes using the built-in class function
%
% Authors: Adam Dai
% Created: 14 Apr 2021
% Updated: 
clear ; clc ;
%% basic 2D 2-ellipsotope
p = 2;
c = zeros(2,1);
G = [2 1; 1 2];

figure(); axis equal; grid on
E = ellipsotope(p,c,G);
plot(E); title('Basic 2D 2-ellipsotope');

%% constrained 2D 2-ellipsotope
p = 2;
c = zeros(2,1);
G = [1 -1 0.2; 1 2 1];
A = [1 -1 1];
b = 0.5;

figure(); axis equal; grid on
E = ellipsotope(p,c,G,A,b);
plot(E); title('Constrained 2D 2-ellipsotope');

%% indexed 2D 2-ellipsotope
p = 2;
c = zeros(2,1);
G = [1 -1 0.2; 1 2 1];
I = {[1 2],3};

figure(); axis equal; grid on
E = ellipsotope(p,c,G,[],[],I);
plot(E); title('Indexed 2D 2-ellipsotope');

%% zonotope using indexed ellipsotope (eqn 8)
p = 2;
c = zeros(2,1);
G = [1 -1 0.2; 1 2 1];
I = {1,2,3};

figure(); axis equal; grid on
E = ellipsotope(p,c,G,[],[],I);
plot(E); title('Zonotope using indexed ellipsotope (eqn 8)')

%% random 2D general 2-ellipsotope
p = 2;
c = zeros(2,1);
G = 2*rand(2,4) - 1;
A = [-1 1 -1 1];
b = 0.5;
I = {[1,2],[3,4]};

figure(); axis equal; grid on
E = ellipsotope(p,c,G,A,b,I);
plot(E); title('Random 2D general ellipsotope');

%% basic 2D 3-ellipsotope
p = 3;
c = zeros(2,1);
G = [2 1; 1 2];

figure(); axis equal; grid on
E = ellipsotope(p,c,G);
plot(E); title('Basic 2D 3-ellipsotope');

%% basic 3D 2-ellipsotope
p = 2;
c = zeros(3,1);
G = eye(3);

figure(); axis equal; grid on; 
E = ellipsotope(p,c,G);
plot(E); title('Basic 3D 2-ellipsotope');

%% general 3D 2-ellipsotope
p = 2;
c = zeros(3,1);
G = [1 0 0 0.5;
     0 1 0 0;
     0 0 1 0.5];
A = [-1 1 -1 0.5];
b = 0.5;
I = {[1,2],[3 4]};

figure(); axis equal; grid on; shading flat
E = ellipsotope(p,c,G,A,b,I);
plot(E, 'facealpha', 1); shading interp; title('General 3D 2-ellipsotope');