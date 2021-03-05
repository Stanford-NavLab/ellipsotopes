%% user parameters
p = 4 ;
% r = 0.5 ;

%% automated from here
% X = make_unit_superellipse_3D(p_norm) ;
X = make_unit_superellipse_2D(p_norm,200) ;

%% plotting=
figure(1) ; clf ; axis equal ; hold on ; grid on ; % view(3)

plot_path(X,'b.')
% patch('faces',F,'vertices',V,'facealpha',0.1)