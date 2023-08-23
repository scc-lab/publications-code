function z_dot = faultyModel(t,z)
    Kp = 4*[1,1,1,1,1]; % x,y,z,phi,theta % Kp = 5, Ki = 2; Kd = 8 for original data set
    Ki = 3*[1,1,1,1,1];
    Kd = 4*[1,1,1,1,1]; % minor anomaly: 8, 4, 2; major anomaly 1: 15,12,2;
    P.Kp = diag(Kp);
    P.Ki = diag(Ki);
    P.Kd = diag(Kd);
    P.xd = [0;0;0];
    % P.tau_phi = 0.00005;
    % P.tau_theta = 0.00005;
    P.tau_psi = 0.00005; % disregard heading angle
    x = z(1:12); % extract x from z
    Ie = z(13:end); % Integral of error
    [u,Ie_dot] = quadControl(P,t,x,Ie); % calculate control input
    x_dot = quadOpenLoop(P,t,x,u); % system model
    z_dot = [x_dot;Ie_dot]; % build z_dot
end