function [A, B, C, K] = obtain_system_matrices( k, trajectory, reach, params )

xnom = trajectory.x_nom; unom = trajectory.u_nom; Knom = trajectory.K_nom;

A = eye(3);
A(1,3) = -unom(1,k-1)*params.dt*sin(xnom(3,k-1));
A(2,3) = unom(1,k-1)*params.dt*cos(xnom(3,k-1));

B = zeros(3,2);
B(1,1) = params.dt*cos(xnom(3,k-1));
B(2,1) = params.dt*sin(xnom(3,k-1));
B(3,2) = params.dt;

K = Knom(:,:,k-1);

C = nan(5,3);
for ibeacon = 1:4
    beacon_pos = reach.beacon_positions(:,ibeacon);
    dnom = sqrt( (xnom(1,k)-beacon_pos(1))^2 + (xnom(2,k)-beacon_pos(2))^2 );
    C(ibeacon,:) = [(xnom(1,k)-beacon_pos(1))/dnom (xnom(2,k)-beacon_pos(2))/dnom 0];
end
C(5,:) = [0 0 1];