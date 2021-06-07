function E_rot = rotated_body(robot, theta, d_h)
    
    w = robot.width; l = robot.length;

    % halfplanes
    z = (w/2)*sin(pi/2-d_h) + (l/2)*sin(d_h);
    x1 = z*[cos(pi/2+theta);sin(pi/2+theta)];
    h1 = [cos(theta+pi/2),sin(theta+pi/2)]; f1 = h1*x1;
    x2 = -z*[cos(pi/2+theta);sin(pi/2+theta)];
    h2 = -[cos(theta+pi/2),sin(theta+pi/2)]; f2 = h2*x2;

    % intersect
    E_rot = halfspace_intersect(robot.circ,h1,f1);
    E_rot = halfspace_intersect(E_rot,h2,f2);

end