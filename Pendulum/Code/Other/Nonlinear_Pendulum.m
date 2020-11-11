clear;

syms omega_n
gValue = 9.81;
mValue = 2;
lValue = 1;
omega_nValue = sqrt(gValue/lValue);
T = 2*pi/omega_nValue;
dtValue = 0.02; 
%%
syms theta(t) theta_d_t(t) th ddth omega_n weight force_th m g l
eqs = [diff(theta)   == theta_d_t;
       diff(theta_d_t) == -omega_n^2*sin(theta)];
   
eqs  = subs(eqs,omega_n,omega_nValue);
vars = [theta, theta_d_t];

[M,F] = massMatrixForm(eqs,vars);

f = M\F;
f_fnc = odeFunction(f, vars);

th0 = [0.1*pi; 0];

tInit  = 0;
tFinal = 40;
ts = [tInit:dtValue:tFinal];

[tout,yout] = ode45(f_fnc,ts,th0);

ddth = -g/l*sin(th);
matlabFunction(ddth,'File','Angular_Accelaration','Vars',[g l th],'Outputs',{'ddth'});
sols_ddth = Angular_Accelaration(gValue, lValue, yout(:,1));

weight = m*g;
weightPlot = double(subs(weight,[m g],[mValue gValue]));
weightVals = double(ones(size(tout)) * weightPlot);

force_th = -weight*sin(th);
matlabFunction(force_th,'File','Force_azimuthal','Vars',[g m th],'Outputs',{'Fth'});
sols_Fth = Force_azimuthal(gValue, mValue, yout(:,1));

%%

noisy_theta = awgn(yout(:,1), 10);
noisy_dtheta = awgn(yout(:,2), 10);
noisy_ddtheta = awgn(sols_ddth, 10);

noisy_weight = awgn(weightVals, 10);
noisy_Fth = awgn(sols_Fth, 10);

figure(1),clf;
hold on

plot(tout,yout(:,1), '-', 'LineWidth', 2);
plot(tout,yout(:,2), '-', 'LineWidth', 2);
plot(tout,sols_ddth, '-', 'LineWidth', 2);

plot(tout, noisy_theta,'.');
plot(tout, noisy_dtheta, '.');
plot(tout, noisy_ddtheta, '.');

title('Non-Linear Motion with Measurements');
grid on;
xlabel('t(s)');
legend('th','dth','ddth','noisy_t_h','noisy_d_t_h','noisy_d_d_t_h')
hold off

figure(2),clf
hold on 
plot(tout,weightVals,'LineWidth',2)
plot(tout,sols_Fth,'LineWidth',2)

plot(tout,noisy_weight,'.')
plot(tout,noisy_Fth,'.')

title('Pendulum Forces')
grid on
xlabel('t(s)')
legend('W','Fth')
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
syms th dth ddth W Fth
syms m g l dt
states = [th; dth; ddth; W; Fth];
vars = [m g l dt th dth ddth W Fth];
equ_th = th + dth*dt + 0.5*dt^2;
equ_dth = dth + ddth*dt;
equ_ddth = -g/l*sin(th);
equ_W = m*g;
equ_Fth = -W*sin(th); 


equ_f = [equ_th; equ_dth; equ_ddth; equ_W; equ_Fth];
equ_F = jacobian(equ_f,states);

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


z = [yout(:,1)'; yout(:,2)'; sols_ddth'; weightVals'; sols_Fth'];
z_noisy = [noisy_theta(2:end)'; noisy_dtheta(2:end)'; noisy_ddtheta(2:end)'; noisy_weight(2:end)'; noisy_Fth(2:end)'];

equ_h = [th; dth; ddth; W; Fth];
equ_H = jacobian(equ_h, states);

matlabFunction(equ_h, 'file','predict_h','Vars', vars);
matlabFunction(equ_H, 'file','predict_hj','Vars', vars);

%% Covariance Matrices

syms cov_th_th  cov_th_dth cov_th_ddth   cov_th_W   cov_th_Fth...
     cov_dth_th cov_dth_dth cov_dth_ddth  cov_dth_W  cov_dth_Fth...
     cov_ddth_th cov_ddth_dth cov_ddth_ddth cov_ddth_W cov_ddth_Fth...
     cov_W_th    cov_W_dth    cov_W_ddth    cov_W_W    cov_W_Fth...
    cov_Fth_th  cov_Fth_dth  cov_Fth_ddth  cov_Fth_W  cov_Fth_Fth

Q = [cov_th_th,   cov_th_dth,   cov_th_ddth,   cov_th_W,   cov_th_Fth;
     cov_dth_th,  cov_dth_dth,  cov_dth_ddth,  cov_dth_W,  cov_dth_Fth;
     cov_ddth_th, cov_ddth_dth, cov_ddth_ddth, cov_ddth_W, cov_ddth_Fth;
     cov_W_th,    cov_W_dth,    cov_W_ddth,    cov_W_W,    cov_W_Fth;
     cov_Fth_th,  cov_Fth_dth,  cov_Fth_ddth,  cov_Fth_W,  cov_Fth_Fth];

Q = double(subs(Q,[cov_th_th,   cov_th_dth,   cov_th_ddth,   cov_th_W,   cov_th_Fth;...
                   cov_dth_th,  cov_dth_dth,  cov_dth_ddth,  cov_dth_W,  cov_dth_Fth;...
                   cov_ddth_th, cov_ddth_dth, cov_ddth_ddth, cov_ddth_W, cov_ddth_Fth;...
                   cov_W_th,    cov_W_dth,    cov_W_ddth,    cov_W_W,    cov_W_Fth;...
                   cov_Fth_th,  cov_Fth_dth,  cov_Fth_ddth,  cov_Fth_W,  cov_Fth_Fth],...
                   [1 0 0 0 0; ...
                    0 1 0 0 0; ...
                    0 0 1 0 0; ...
                    0 0 0 1 0; ...
                    0 0 0 0 1]));


R = [cov_th_th,   cov_th_dth,   cov_th_ddth,   cov_th_W,   cov_th_Fth;
     cov_dth_th,  cov_dth_dth,  cov_dth_ddth,  cov_dth_W,  cov_dth_Fth;
     cov_ddth_th, cov_ddth_dth, cov_ddth_ddth, cov_ddth_W, cov_ddth_Fth;
     cov_W_th,    cov_W_dth,    cov_W_ddth,    cov_W_W,    cov_W_Fth;
     cov_Fth_th,  cov_Fth_dth,  cov_Fth_ddth,  cov_Fth_W,  cov_Fth_Fth];
 
R = double(subs(R,[cov_th_th,   cov_th_dth,   cov_th_ddth,   cov_th_W,   cov_th_Fth;...
                   cov_dth_th,  cov_dth_dth,  cov_dth_ddth,  cov_dth_W,  cov_dth_Fth;...
                   cov_ddth_th, cov_ddth_dth, cov_ddth_ddth, cov_ddth_W, cov_ddth_Fth;...
                   cov_W_th,    cov_W_dth,    cov_W_ddth,    cov_W_W,    cov_W_Fth;...
                   cov_Fth_th,  cov_Fth_dth,  cov_Fth_ddth,  cov_Fth_W,  cov_Fth_Fth],...
                   [1 0 0 0 0; ...
                    0 1 0 0 0; ...
                    0 0 1 0 0; ...
                    0 0 0 1 0; ...
                    0 0 0 0 1]));

%%
% state covariance matrix P
P_int = [pi/180 0 0 0 0;...
         0 pi/180 0 0 0;...
         0 0 pi/180 0 0;...
         0 0 0 1 0;...
         0 0 0 0 1];
x_int = [0.1*pi; 0; -9.81*0.1*pi; 2*9.81; -5];

% initialise

P = P_int;
x = x_int;
x_predict = x_int;
xs = zeros(5,length(ts));
xh = zeros(5,length(ts));
I = eye(5);
x_hat = zeros(5,1);


%% Implement the Kalman Filter

dts = 1:1:length(ts)-1; %start indexing at t1 for measurements z
%%
xs(:,1) = x_int;
xh(:,1) = x_int;
Ph = [diag(P_int)]; % line 233

tstart = tic;
for i = dts
    % predict
    x_hat = predict_f(mValue,gValue,lValue,dtValue,x(1),x(2),x(3),x(4),x(5));
    F = predict_fj(mValue,gValue,lValue,dtValue,x(1),x(2),x(3),x(4),x(5));
    P_hat = F*P*F' + Q;
    
    % update
    H = predict_hj(mValue,gValue,lValue,dtValue,x(1),x(2),x(3),x(4),x(5));
    S = H*P_hat+H' + R;
    K = P_hat*H'/S;
    y = z_noisy(:,i) - H*x_hat;
    x = x_hat + K*y;
    P = (I-K*H)*P_hat;

    xs(:,i+1) = x;
    
    Ph = [Ph, diag(P)]; 
    x_predict = predict_f(mValue,gValue,lValue,dtValue,x_predict(1),x_predict(2),x_predict(3),x_predict(4),x_predict(5));
    xh(:,i+1) = x_predict;
end
telapsed = toc(tstart);

%%
% plot xh (predicted output of states)
figure(3),clf
plot(ts, xh);
legend(["th" "dth" "ddth" "W" "Fth"]);
title('Predicted Output of States')
grid on
xlabel('t(s)')

% plot xs (Kalman Estimate)
figure(4),clf
plot(ts, xs);
legend(["th" "dth" "ddth" "W" "Fth"]);
title('Kalman Filter Output')
grid on
xlabel('t(s)')


xs = xs';
z = z';

% plot true state vs kalman estimate
figure(5),clf
subplot(2,2,1)
hold on
plot(ts, xs(:,1));
plot(ts, z(:,1),'LineWidth',1.5);
legend(["th", "th-true"])
title('Position')
grid on
xlabel('t(s)')

subplot(2,2,2)
hold on
plot(ts, xs(:,2));
plot(ts, z(:,2),'LineWidth',1.5);
legend(["dth", "dth-true"])
title('Velocity')
grid on
xlabel('t(s)')

subplot(2,2,3)
hold on
plot(ts, xs(:,3));
plot(ts, z(:,3),'LineWidth',1.5);
title('Acceleration')
legend(["ddth", "ddth-true"])
grid on
xlabel('t(s)')

subplot(2,2,4)
hold on
plot(ts, xs(:,5));
plot(ts, z(:,5),'LineWidth',1.5);
title('Fth')
legend(["F", "F-true"])
grid on
xlabel('t(s)')


%%
RMS  = [(xs(:,1) - z(:,1)).^2, ...
		(xs(:,2) - z(:,2)).^2, ...
		(xs(:,3) - z(:,3)).^2, ...
		(xs(:,4) - z(:,4)).^2, ...
		(xs(:,5) - z(:,5)).^2];

RMS_ave = sum(RMS);
RMS_ave = RMS_ave./length(ts)

% uncorellated estimate error should similar to RMS values

%%

figure(6),clf
subplot(2,2,1)
plot(ts, Ph(1,:));
legend("P-th")


subplot(2,2,2)
plot(ts, Ph(2,:));
legend("P-dth")


subplot(2,2,3)
plot(ts, Ph(3,:));
legend("P-ddth")


subplot(2,2,4)
plot(ts, Ph(5,:));
legend("P-F")

