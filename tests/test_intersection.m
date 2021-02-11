% Intersect 2 circles represented as ellipsotopes

E1 = ellipsotope(2,[0;0],eye(2));
E2 = ellipsotope(2,[1;0],eye(2));

E_int = E1 & E2;


figure(1) ; clf ; axis equal ; hold on ; grid on ;

plot(E1);
plot(E2);
plot(E_int);