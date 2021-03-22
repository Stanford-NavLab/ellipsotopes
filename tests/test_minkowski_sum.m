% Minkowski sum of 2 circles

E1 = ellipsotope(2,[0;0],diag([2,1]));
E2 = ellipsotope(2,[1;0],diag([1,2]));

%E_sum = E1 + E2;
E_sum = ellipsotope(2,[1;0],[diag([2,1]) diag([1,2])]);

figure(1) ; clf ; axis equal ; hold on ; grid on ;

plot(E_sum);
plot(E1);
plot(E2);
