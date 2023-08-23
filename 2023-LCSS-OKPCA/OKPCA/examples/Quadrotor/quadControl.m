function [u,Ie_dot] = quadControl(P,t,x,Ie)
    phi = x(7);
    theta = x(8);
    g = 9.81;
    m = 0.579902;
    % Command angles
    ux = x(4);
    uy = x(5);
    uz = x(6);
    F = m*(g-uz)/(cos(phi)*cos(theta));
    phi_com = atan(uy*cos(theta)/(g-uz));
    theta_com = atan(ux/(uz-g));
    e = [P.xd - x(1:3); phi_com - phi; theta_com - theta]; % x, y, z, phi, theta
    e_dot = [0 - x(4:6); 0-x(10:11)]; % u, v, w, phi_dot, theta_dot
    Ie_dot = e; 
    u = P.Kp*e + P.Ki*Ie + P.Kd*e_dot;
end