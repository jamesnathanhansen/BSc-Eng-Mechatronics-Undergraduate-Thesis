%% Simple Pendulum EOM 
% Simulation of Pendulum Motion
% Create the angular accelaration EOM
syms m a g l theta(t) omega_n dt
eqn = m*a == -m*g*sin(theta);
eqn = subs(eqn,a,l*diff(theta,2));
eqn = isolate(eqn,diff(theta,2));
eqn = subs(eqn,g/l,omega_n^2);
assume(omega_n,'real')

% Linearize
syms x
eqnLinear = subs(eqn,sin(theta(t)),theta(t));

syms theta_0 dtheta_0
dtheta = diff(theta);
cond = [theta(0) == theta_0, dtheta(0) == dtheta_0];

% Solve analytically for angular position, velocity, acceleration
thetaSolLin(t) = dsolve(eqnLinear,cond);
dthetaSolLin(t) = diff(thetaSolLin(t));
ddthetaSolLin(t) = -omega_n^2*thetaSolLin;

% Parameters
gValue = 9.81;
lValue = 1;
mValue = 1;
omega_nValue = sqrt(gValue/lValue);
dtValue = 0.02; 

% Inital Conditions
theta_0Value  = deg2rad(10); % Solution only valid for small angles.
dtheta_0Value = 0; % Initially at rest.
ddtheta_0Value = 0;

% Sub in Parameters
vars   = [omega_n      theta_0      dtheta_0];
values = [omega_nValue theta_0Value dtheta_0Value];
thetaSolPlot = subs(thetaSolLin,vars,values);

vars   = [omega_n      theta_0      dtheta_0];
values = [omega_nValue theta_0Value dtheta_0Value];
dthetaSolPlot = subs(dthetaSolLin,vars,values);

vars   = [omega_n      theta_0      dtheta_0];
values = [omega_nValue theta_0Value dtheta_0Value];
ddthetaSolPlot = subs(ddthetaSolLin,vars,values);

% Create vector of output for time ts
ts = 0:dtValue:25;

% Find values for angular position, velocity and acceleration
thetaSolVal = double(thetaSolPlot(ts));
dthetaSolVal = double(dthetaSolPlot(ts));
ddthetaSolVal = double(ddthetaSolPlot(ts));

motion = [thetaSolVal; dthetaSolVal; ddthetaSolVal];

% Add Sensor noise
reset(RandStream.getGlobalStream) % Produce reproducible results each run
noisy_thetaSolVal = awgn(thetaSolVal,10,'measured');
noisy_dthetaSolVal = awgn(dthetaSolVal,10,'measured');
noisy_ddthetaSolVal = awgn(ddthetaSolVal,10,'measured');

motion_noise = [noisy_thetaSolVal; noisy_dthetaSolVal; noisy_ddthetaSolVal];
%% Plotting
figure(1),clf
grid on 
title('Linear Pendulum Motion');
xlabel('t (s)');
hold on

% Plot EOM
plot(ts,thetaSolVal,'-','LineWidth',1)
plot(ts,dthetaSolVal,'-','LineWidth',1)
plot(ts,ddthetaSolVal,'-','LineWidth',1)

% Plot with noise
plot(ts,noisy_thetaSolVal,'.')
plot(ts,noisy_dthetaSolVal,'.')
plot(ts,noisy_ddthetaSolVal,'.')

legend('Angular Position','Angular Velocity','Angular Accelaration','Angular Position with Noise','Angular Velocity with Noise','Angular Accelaration with Noise')
hold off
%% Kalman Filter Model
%% Prediction Step Design
% State Transisition Function
F = [1,   dt, 0.5*dt^2;...
     0,   1,  dt;...
    -g/l, 0,  0];

F = double(subs(F,[dt g l], [dtValue gValue lValue]));

% Process Noise Matrix
Q = [0 0 0; 
     0 0 0; 
     0 0 0]; 
%% Update Step Design
% Measurements
z = [noisy_thetaSolVal(2:end); noisy_dthetaSolVal(2:end); noisy_ddthetaSolVal(2:end)];

% Measuremnt Matrix
H = [1, 0, 0; ...
     0, 1, 0;...
     0, 0, 1];

% Find Measurement Noise Variances for States
diff_theta_sq = ((thetaSolVal - noisy_thetaSolVal).^2);
sum_theta = sum(diff_theta_sq);
variance_theta = sum_theta/length(thetaSolVal);

diff_dtheta_sq = ((dthetaSolVal - noisy_dthetaSolVal).^2);
sum_dtheta = sum(diff_dtheta_sq);
variance_dtheta = sum_dtheta/length(dthetaSolVal);

diff_ddtheta_sq = ((ddthetaSolVal - noisy_ddthetaSolVal).^2);
sum_ddtheta = sum(diff_ddtheta_sq);
variance_ddtheta = sum_ddtheta/length(ddthetaSolVal);

% Measurement Noise Matrix
R = [variance_theta, 0,               0;
     0,              variance_dtheta, 0;
     0,              0,               variance_ddtheta];
%% Initialise KF
% State Covariance Matrix P
P_int = R.*10;
x_int = [deg2rad(14); 0; 0];

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
    x_hat = F*x;
    P_hat = F*P*F' + Q;
    
    % Update Step
    S = H*P_hat+H' + R;
    K = P_hat*H'/S;
    y = z(:,i) - H*x_hat;
    x = x_hat + K*y;
    P = (I-K*H)*P_hat;
    
    % Append Current Values to Output Vector
    xs(:,i+1) = x;
    P_cov = [P_cov, diag(P)];
    
    % Prediction Step w/o Updates
    x_predict = F*x_predict;
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
title(tiled, 'True States vs KF Estimated States vs Measurements')
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