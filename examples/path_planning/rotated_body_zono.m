function Z_bound = rotated_body_zono(robot, theta, d_h)
    
    w = robot.width; l = robot.length;
    r = norm([w l]) / 2;
    
    z = (w/2)*sin(pi/2-d_h) + (l/2)*sin(d_h);
    psi = atan2(w,l);
    if d_h < psi
        H = r * cos(psi - d_h);
    else
        H = r;
    end

    % creating bounding zonotope
    G_bound = diag([H,z]);
    Z_bound = rotation_matrix_2D(theta) * zonotope([0;0],G_bound);

end