clear;
syms omega_n
gValue = 9.81;
mValue = 2;
lValue = 1;
omega_nValue = sqrt(gValue/lValue);
T = 2*pi/omega_nValue;
dtValue = 0.02; 
%%
syms theta(t) theta_d_t(t) th ddth omega_n  m g l
eqs = [diff(theta)   == theta_d_t;
       diff(theta_d_t) == -omega_n^2*sin(theta)];
   
eqs  = subs(eqs,omega_n,omega_nValue);
vars = [theta, theta_d_t];

[M,F] = massMatrixForm(eqs,vars);

f = M\F;
f_fnc = odeFunction(f, vars);

th0 = [0.1*pi; 0];

tInit  = 0;
tFinal = 50;
ts = [tInit:dtValue:tFinal];

[tout,yout] = ode45(f_fnc,ts,th0);

ddth = -g/l*sin(th);
matlabFunction(ddth,'File','Angular_Accelaration','Vars',[g l th],'Outputs',{'ddth'});
sols_ddth = Angular_Accelaration(gValue, lValue, yout(:,1));

%%
figure(1),clf;
hold on

plot(tout,yout(:,1), '-', 'LineWidth', 2);
plot(tout,yout(:,2), '-', 'LineWidth', 2);
plot(tout,sols_ddth, '-', 'LineWidth', 2);

noisy_theta = awgn(yout(:,1), 10);
noisy_dtheta = awgn(yout(:,2), 10);
noisy_ddtheta = awgn(sols_ddth, 10);

plot(tout, noisy_theta,'.');
plot(tout, noisy_dtheta, '.');
plot(tout, noisy_ddtheta, '.');

title('Non-Linear Motion with Measurements');
grid on;
xlabel('t(s)');
legend('th','dth','ddth','noisy_t_h','noisy_d_t_h','noisy_d_d_t_h')
hold off
%% Extended Kalman Filter Model
% Design the state variables
% We want to track the mass of the pendulum assuming the string is light and inextensible and no friction.
% This means we need state variables that will tell us about the mass.
% So we need weight (W) (knowing the weight means knows the mass). 

% equ_W = m*g;

% Similarly, we need a way to describe W in terms of something measureble.
% We know that force due to gravity results in reaction force in the radial
% direction Fr and the tangential direction Fth along the strong. 

% equ_Fth = -W*sin(th);

% What tells about th? dth.

% Let's assume we have a gyroscope. This can measure angular velocity. We
% can integrate this to get th.

% equ_th = th0 + dth*dt

% Fth and sin(th). So Fth and th are states as well. 
%  What tells us about dth? ddth. So dth and ddth
% are states as well. 

% states
% x = [th; dth; ddth; W; Fth]

% x_hat = f(x,u) = [equ_th; equ_dth; equ_ddth; equ_W; equ_Fth]

% equ_th = th0 + dth*dt
% equ_dth = dth0 + ddth*dt
% equ_ddth = ddth0 - g/l*sin(th)
% equ_W = m*g;
% equ_Fth = -W*sin(th);

%% Prediction
syms th dth ddth
syms m g l dt
states = [th; dth; ddth];
vars = [m g l dt th dth ddth];
equ_th = th + dth*dt;
equ_dth = dth + ddth*dt;
equ_ddth = - g/l*sin(th);

equ_f = [equ_th; equ_dth; equ_ddth]
equ_F = jacobian(equ_f,states)

matlabFunction(equ_f, 'file','predict_f','Vars', vars);
matlabFunction(equ_F, 'file','predict_fj','Vars', vars);
%% Measuremnts

% In terms of measurements, we have to ask what equipment do we have that
% can measure states in the real world as well? 
% Do we have enough states to do the update step? 
% Assuming we have direct access to to any state, what measurements are
% required to fully describe my system?  


% z = [th; dth; ddth; W; Fth]

% h(x,u) needs to translate the predicted states to measurements
% [5x1] = [5x1] - [5x5][5x1]
% y = z - h(x,u)
% h(x,u) = [th; dth; ddth; W; Fth]
% H = [1 0 0 0 0; 
%      0 1 0 0 0;
%      0 0 1 0 0;
%      0 0 0 1 0;
%      0 0 0 0 1]


z = [yout(:,1)'; yout(:,2)'; sols_ddth'];
z_noisy = [noisy_theta(2:end)'; noisy_dtheta(2:end)'; noisy_ddtheta(2:end)'];

equ_h = [th; dth; ddth]
equ_H = jacobian(equ_h, states)

matlabFunction(equ_h, 'file','predict_h','Vars', vars);
matlabFunction(equ_H, 'file','predict_hj','Vars', vars);

%% Covariance Matrices

syms cov_th_th  cov_th_dth cov_th_ddth...
     cov_dth_th cov_dth_dth cov_dth_ddth...
     cov_ddth_th cov_ddth_dth cov_ddth_ddth...


Q = [cov_th_th,   cov_th_dth,   cov_th_ddth;
     cov_dth_th,  cov_dth_dth,  cov_dth_ddth;
     cov_ddth_th, cov_ddth_dth, cov_ddth_ddth];

Q = double(subs(Q,[cov_th_th,   cov_th_dth,   cov_th_ddth;...
                   cov_dth_th,  cov_dth_dth,  cov_dth_ddth;...
                   cov_ddth_th, cov_ddth_dth, cov_ddth_ddth],...
                   [1 0 0; ...
                    0 1 0; ...
                    0 0 1]))
Q = Q./100;

R = [cov_th_th,   cov_th_dth,   cov_th_ddth;
     cov_dth_th,  cov_dth_dth,  cov_dth_ddth;
     cov_ddth_th, cov_ddth_dth, cov_ddth_ddth];

R = double(subs(Q,[cov_th_th,   cov_th_dth,   cov_th_ddth;...
                   cov_dth_th,  cov_dth_dth,  cov_dth_ddth;...
                   cov_ddth_th, cov_ddth_dth, cov_ddth_ddth],...
                   [1 0 0; ...
                    0 1 0; ...
                    0 0 1]))
R = R./100;
%%
% state covariance matrix P
P_int = [pi/180 0 0; 0 pi/180 0; 0 0 pi/180];
x_int = [0.1*pi; 0; -10*0.1*pi];

% initialise

P = P_int;
x = x_int;
x_predict = x_int;
xs = zeros(3,length(ts));
xh = zeros(3,length(ts));
xm = zeros(3,length(ts));
I = eye(3);
x_hat = zeros(3,1);


%% Implement the Kalman Filter

dts = 1:1:length(ts)-1; %start indexing at t1 for measurements z
%%
xs(:,1) = x_int;
xh(:,1) = x_int;
xm(:,1) = x_int;

tstart = tic;
for i = dts
    % predict
    x_hat = predict_f(mValue,gValue,lValue,dtValue,x(1),x(2),x(3));
    F = predict_fj(mValue,gValue,lValue,dtValue,x(1),x(2),x(3));
    P_hat = F*P*F' + Q;
    
    % update
    H = predict_hj(mValue,gValue,lValue,dtValue,x(1),x(2),x(3));
    S = H*P_hat+H' + R;
    K = P_hat*H'/S;
    y = z_noisy(:,i) - H*x_hat;
    x = x_hat + K*y;
    P = (I-K*H)*P_hat;
    P = (P+P')/2;
    
    % kalman estimate
    xs(:,i+1) = x;
    
    % predict output with kalman adjustment
    xh(:,i+1) = x_hat;
    
    % model output without kalman adjustment
    x_predict = predict_f(mValue,gValue,lValue,dtValue,x_predict(1),x_predict(2),x_predict(3));
    xm(:,i+1) = x_predict;
end
telapsed = toc(tstart)

%%
% plot xm 
figure(3),clf
plot(ts, xm);
legend(["th" "dth" "ddth"]);
title('Model Output without Kalman Adjustment')
grid on
xlabel('t(s)')

% plot xh (Predict)
figure(4),clf
plot(ts, xh);
legend(["th" "dth" "ddth"]);
title('Predict Output')
grid on
xlabel('t(s)')

% plot xs (Kalman Estimate)
figure(5),clf
plot(ts, xs);
legend(["th" "dth" "ddth"]);
title('Kalman Filter Output')
grid on
xlabel('t(s)')