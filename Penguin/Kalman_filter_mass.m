function [ts,xs,xp,z,Pvals] = Kalman_filter_mass(Q_cov,R_cov,P_int_cov,x_int_val,measurements)
    %% Design the state variables x
    syms m g

    states = [m];
    vars = [m];
    %% Design the state transition model F
    equ_m = m;
    equ_f = [equ_m];
    equ_F = jacobian(equ_f,states);

    matlabFunction(equ_f, 'file','predict_f','Vars', vars);
    matlabFunction(equ_F, 'file','predict_fj','Vars', vars);
    %% Design the observation model H
    equ_h = [m];
    equ_H = jacobian(equ_h, states);

    matlabFunction(equ_h, 'file','predict_h','Vars', vars);
    matlabFunction(equ_H, 'file','predict_hj','Vars', vars);
    %% Design the state covaraince P
    P = P_int_cov;
    %% Design the process noise covariance Q
    Q = Q_cov;
    %% Design the measurement noise covariance R
    R = R_cov;
    %% Values
    mValue = x_int_val;
    %% Measurements
    z = measurements; % Units in kg
    %% Initialise
    I = eye(1);
    x_int = mValue; % Convert to N
    x = x_int; % Initial Weight
    x_predict = x_int;

    %% Time Parameters
    syms f dt 
    f = 50; % Hz
    dt = 1/f;
    tInt = 0; 
    tFinal = dt*length(z);
    ts = tInt:dt:tFinal;

    %% Vectors for Outputs, Plots etc.

    xs = zeros(1,length(ts));
    xs(:,1) = x_int;

    xp = zeros(1,length(ts));
    xp(:,1) = x_int;

    Pvals = zeros(1,length(ts));
    Pvals(:,1) = P;
    
    %% Run Algorithm
    dts = 1:1:length(ts)-1;

    for i = dts
        % Prediction Step
        x_hat = predict_f(x);
        F = predict_fj(x);
        P_hat = F*P*F' + Q;

        % Update Step
        H = predict_hj(x);
        S = H*P_hat+H' + R;
        K = P_hat*H'/S;
        y = z(i) - H*x_hat;
        x = x_hat + K*y;
        P = (I-K*H)*P_hat;

        % Append Current Values to Output Vectors
        xs(:,i+1) = x;
        Pvals(:,i+1) = P;

        % Prediction Step w/o Updates
        x_predict = predict_f(x_predict);
        xp(:,i+1) = x_predict;
    end
end