% Minkowski sum of 2 circles

% c1 = [0;0]; G1 = [1 0 0 0; 0 1 0 0];
% J1 = {[1,2],[3,4]}; A1 = [1 0 -1 0; 0 1 0 -1]; b1 = [1; 0];
% E1 = ellipsotope(2, c1, G1, J1, A1, b1);

% c1 = [0;0]; G1 = [1 0 0 0; 0 2 0 0];
% J1 = {[1,2],[3,4]}; A1 = [1 0 -1 0; 0 2 0 -1]; b1 = [1; 0];
% E2 = ellipsotope(2, c2, G2, J2, A2, b2);

E_temp1 = ellipsotope(2,[0;0],eye(2));
E_temp2 = ellipsotope(2,[1;0],eye(2));
E1 = E_temp1 & E_temp2;

E_temp1 = ellipsotope(2,[0;0],diag([1,2]));
E_temp2 = ellipsotope(2,[1;0],eye(2));
E2 = E_temp1 & E_temp2;
E2 = E2 + [1;1];

%E_sum = E1 + E2;
%E_sum = ellipsotope(2,[1;0],[diag([2,1]) diag([1,2])]);

figure(1) ; clf ; axis equal ; hold on ; grid on ;

%plot(E_sum);
plot(E1);
plot(E2);
