function x_dot = quadOpenLoop(P,t,x,u)
% Simplified quadrotor model, small-angle approx. & coriolis terms neglected
    g = 9.81;
    m = 0.579902;
    Jx = 0.002261; % kg*m^2
    Jy = 0.002824; % kg*m^2
    Jz = 0.002097; % kg*m^2
    
    % Control inputs
    ux = u(1);
    uy = u(2);
    uz = u(3);
    tau_phi = u(4);
    tau_theta = u(5);
    p_dot = 1/Jx * tau_phi;
    q_dot = 1/Jy * tau_theta;
    r_dot = 1/Jz * P.tau_psi;
    
    x_dot = [x(4); % xdot
          x(5); % ydot
          x(6); % zdot
          ux; % u
          uy; % v
          uz; % w
          x(10); % phidot
          x(11); % thetadot
          x(12); % psidot
          p_dot;
          q_dot;
          r_dot];
end