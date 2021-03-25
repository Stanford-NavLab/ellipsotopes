% Minkowski sum of 2 circles

c1 = [0;0]; G1 = diag([2,1]);
J1 = {1,2}; A1 = [1 0; 0 0.5]; b1 = [2; 2];
E1 = ellipsotope(2, c1, G1, J1, A1, b1);

c2 = [1;0]; G2 = diag([1,2]);
J2 = {1,2}; A2 = [1 1; 0 2]; b2 = [1; -2];
E2 = ellipsotope(2, c2, G2, J2, A2, b2);

%E_sum = E1 + E2;
%E_sum = ellipsotope(2,[1;0],[diag([2,1]) diag([1,2])]);

figure(1) ; clf ; axis equal ; hold on ; grid on ;

%plot(E_sum);
plot(E1);
plot(E2);
