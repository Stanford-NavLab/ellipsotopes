%% description
% Robot path planning example
% Robot's volume is represented by a zonotopic ellipsotope, and the
% (probabilistically-bounded) uncertain position of the robot is
% represented by an ellipsoidal ellipsotope.
%
% Robot uses Kalman filter to estimate its position. KF state estimate and
% covariance is used to create uncertain position ellipse.

%% parameters

dt = 0.1;

% state 
x0 = zeros(4,1);

A = [1 0 dt 0; 
     0 1 0 dt;
     0 0 1 0;
     0 0 0 1];
 
B = [