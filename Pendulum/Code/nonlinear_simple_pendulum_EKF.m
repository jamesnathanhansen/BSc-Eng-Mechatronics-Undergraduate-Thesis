%% Nonlinear Simple Pendulum EOM 
% Simulation of Pendulum Motion
% Create the angular accelaration EOM
clear;

% Parameters
gValue = 9.81;
mValue = 2;
lValue = 1;
omega_nValue = sqrt(gValue/lValue);
dtValue = 0.02;

% Rewrite the second-order ODE as a system of first-order ODEs
syms theta(t) dtheta(t) omega_n
eqs = [diff(theta)   == dtheta;
       diff(dtheta) == -omega_n^2*sin(theta)];
   
eqs  = subs(eqs,omega_n,omega_nValue);
vars = [theta, dtheta];

% Find the mass matrix M of the system & the right sides of the equations F
[M,F] = massMatrixForm(eqs,vars);

% Simplify further computations
f = M\F;

% Convert f to a MATLAB function handle by using odeFunction. The resulting
% function handle is the input to the MATLAB ODE solver ode45.
f_fnc = odeFunction(f, vars);

%  Store the initial conditions of θ and dθ/dt in the variable x0.
th0 = [deg2rad(20); 0];

% Specify a time interval from 0 s to 10 s for finding the solution.
tInit  = 0;
tFinal = 40;
ts = tInit:dtValue:tFinal;

% Solve the ODE for values of angular position and velocity 
[tout,yout] = ode45(f_fnc,ts,th0);

% Create function for angular accelaration
syms th
ddth = -omega_n^2*sin(th);
matlabFunction(ddth,'File','Angular_Accelaration','Vars',[omega_n th],'Outputs',{'ddth'});

% Solve for values of angular acceleration
sols_ddth = Angular_Accelaration(omega_nValue, yout(:,1));

motion = [yout(:,1)';yout(:,2)';sols_ddth'];

% Add Sensor noise
reset(RandStream.getGlobalStream) % Produce reproducible results each run
noisy_theta = awgn(yout(:,1),10,'measured');
noisy_dtheta = awgn(yout(:,2),10,'measured');
noisy_ddtheta = awgn(sols_ddth,10,'measured');

motion_noise = [noisy_theta';noisy_dtheta';noisy_ddtheta'];
%% Plotting
% Plot EOM
figure(1),clf;
grid on
title('Nonlinear Pendulum Motion');
xlabel('t (s)');
hold on

plot(tout,yout(:,1),'LineWidth',1);
plot(tout,yout(:,2),'LineWidth',1);
plot(tout,sols_ddth,'LineWidth',1);

% Plot EOM with noise
plot(tout, noisy_theta,'.');
plot(tout, noisy_dtheta, '.');
plot(tout, noisy_ddtheta, '.');

legend('Angular Position','Angular Velocity','Angular Accelaration','Angular Position with Noise','Angular Velocity with Noise','Angular Accelaration with Noise')
hold off
%% Extended Kalman Filter Model
%% Prediction Step Design
% State Transistion Function
syms th dth ddth
syms m g l dt

% States
states = [th; dth; ddth];
vars = [m g l dt th dth ddth];

% State Transistion Equations f
equ_th = th + dth*dt + 0.5*dt^2;
equ_dth = dth + ddth*dt;
equ_ddth = -g/l*sin(th);

equ_f = [equ_th; equ_dth; equ_ddth];

% Discrete State Transition Matrix F
equ_F = jacobian(equ_f,states);

% Write to MATLAB function 
matlabFunction(equ_f, 'file','predict_f','Vars', vars);
matlabFunction(equ_F, 'file','predict_fj','Vars', vars);

% Process Noise Matrix
Q = [1 0 0; 
     0 1 0; 
     0 0 1];
 
Q = Q./1000;
%% Update Step Design
% Measurements
z = [noisy_theta(2:end)'; noisy_dtheta(2:end)'; noisy_ddtheta(2:end)'];

% Measurement Equations h
equ_h = [th; dth; ddth];

% Measurement Matrix H
equ_H = jacobian(equ_h, states);

% Write to MATLAB function
matlabFunction(equ_h, 'file','predict_h','Vars', vars);
matlabFunction(equ_H, 'file','predict_hj','Vars', vars);

% Find Measurement Noise Variances for States
diff_theta_sq = ((yout(:,1) - noisy_theta).^2);
sum_theta = sum(diff_theta_sq);
variance_theta = sum_theta/length(yout);

diff_dtheta_sq = ((yout(:,2) - noisy_dtheta).^2);
sum_dtheta = sum(diff_dtheta_sq);
variance_dtheta = sum_dtheta/length(yout(:,2));

diff_ddtheta_sq = ((sols_ddth - noisy_ddtheta).^2);
sum_ddtheta = sum(diff_ddtheta_sq);
variance_ddtheta = sum_ddtheta/length(sols_ddth);

% Measurement Noise Matrix
R = [variance_theta, 0,               0;
     0,              variance_dtheta, 0;
     0,              0,               variance_ddtheta];
%% Initialise KF
% State covariance matrix P
P_int = R.*10;
x_int = [deg2rad(89); 0; 0];

P = P_int;
x = x_int;
x_predict = x_int;

xs = zeros(3,length(ts));
xp = zeros(3,length(ts));
x_hat = zeros(3,1);

xs(:,1) = x_int;
xp(:,1) = x_int;
P_cov = [diag(P_int)];

I = eye(3);
%% Implement the Kalman Filter
dts = 1:1:length(ts)-1;

for i = dts
    % Prediction Step
    F = predict_fj(mValue,gValue,lValue,dtValue,x(1),x(2),x(3));
    x_hat = predict_f(mValue,gValue,lValue,dtValue,x(1),x(2),x(3));
    P_hat = F*P*F' + Q;
    
    % Update Step
    H = predict_hj(mValue,gValue,lValue,dtValue,x(1),x(2),x(3));
    S = H*P_hat+H' + R;
    K = P_hat*H'/S;
    y = z(:,i) - H*x_hat;
    x = x_hat + K*y;
    P = (I-K*H)*P_hat;
    
    % Append Current Values to Output Vector
    xs(:,i+1) = x;
    P_cov = [P_cov, diag(P)];
    
    % Prediction Step w/o Updates 
    x_predict = predict_f(mValue,gValue,lValue,dtValue,x_predict(1),x_predict(2),x_predict(3));
    xp(:,i+1) = x_predict;
end

%% Plotting
% Plot xp (Predicted Output of States w/o Updates)
figure(2),clf
plot(ts, xp,'LineWidth',1);
grid on
legend(["th" "dth" "ddth"]);
title('Predicted Output of States')
xlabel('t (s)')

% Plot xs (Kalman Estimate)
figure(3),clf
plot(ts, xs,'LineWidth',1);
grid on
legend(["th" "dth" "ddth"]);
title('Kalman Filter Output')
xlabel('t (s)')

% Plot P_cov (Error Covariances)
figure(4),clf
tiled = tiledlayout('flow');
title(tiled, 'State Error Covaraince')
xlabel(tiled, 't (s)')

nexttile
plot(ts, P_cov(1,:),'LineWidth',1);
grid on
title("Variance in Position ")

nexttile
plot(ts, P_cov(2,:),'LineWidth',1);
grid on
title("Variance in Velocity")

nexttile
plot(ts, P_cov(3,:),'LineWidth',1);
grid on
title("Variance in Accelaration")

% Plot Predicted States
figure(5),clf
tiled = tiledlayout('flow');
title(tiled, 'Predicted Output of States without KF Updates vs True States')
xlabel(tiled, 't (s)')

nexttile
hold on
plot(ts, motion(1,:),'LineWidth',1);
plot(ts, xp(1,:),'LineWidth',1);
hold off
grid on
legend(["True th" "Predicted th"]);
title('Angular Position')
ylabel('${\theta}$ (rad)','interpreter','latex')

nexttile
hold on
plot(ts, motion(2,:),'LineWidth',1);
plot(ts, xp(2,:),'LineWidth',1);
hold off
grid on
legend(["True dth" "Predicted dth"]);
title('Angular Velocity')
ylabel('$\dot{\theta}$ (rad/s)','interpreter','latex')

nexttile
hold on
plot(ts, motion(3,:),'LineWidth',1);
plot(ts, xp(3,:),'LineWidth',1);
hold off
grid on
legend(["True ddth" "Predicted ddth"]);
title('Angular Accelaration')
ylabel('$\ddot{\theta}$ ${(rad/s^2)}$','interpreter','latex')

% Plot KF Estimated States
figure(6),clf
tiled = tiledlayout('flow');
title(tiled, 'True States vs KF Estimated States vs Measurements [Q=1]')
xlabel(tiled, 't (s)')

nexttile
hold on
plot(ts, motion(1,:),'LineWidth',1);
plot(ts, xs(1,:),'LineWidth',1);
plot(ts, motion_noise(1,:),'.');
hold off
grid on
legend(["th" "KF th" "measured th"]);
title('Angular Position')
ylabel('${\theta}$ (rad)','interpreter','latex')

nexttile
hold on
plot(ts, motion(2,:),'LineWidth',1);
plot(ts, xs(2,:),'LineWidth',1);
plot(ts, motion_noise(2,:),'.');
hold off
grid on
legend(["dth" "KF dth" "measured dth"]);
title('Angular Velocity')
ylabel('$\dot{\theta}$ (rad/s)','interpreter','latex')

nexttile
hold on
plot(ts, motion(3,:),'LineWidth',1);
plot(ts, xs(3,:),'LineWidth',1);
plot(ts, motion_noise(3,:),'.');
hold off
grid on
legend(["ddth" "KF ddth" "measured ddth"]);
title('Angular Accelaration')
ylabel('$\ddot{\theta}$ ${(rad/s^2)}$','interpreter','latex')