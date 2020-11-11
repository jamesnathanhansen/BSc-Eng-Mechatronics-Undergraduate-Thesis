% Parameters
gValue = 9.81;
rValue = 1;
omega_0Value = sqrt(gValue/rValue);

% Rewrite the second-order ODE as a system of first-order ODEs.
syms theta(t) theta_t(t) omega_0
eqs = [diff(theta)   == theta_t;
       diff(theta_t) == -omega_0^2*sin(theta)]
   
eqs  = subs(eqs,omega_0,omega_0Value);
vars = [theta, theta_t];

% Find the mass matrix M of the system and the right sides of the equations
% F.
[M,F] = massMatrixForm(eqs,vars)

% Simplify further computations
f = M\F

% Convert f to a MATLAB function handle by using odeFunction. The resulting
% function handle is the input to the MATLAB ODE solver ode45.
f = odeFunction(f, vars)

%  Store the initial conditions of θ and dθ/dt in the variable x0.
x0 = [0.1*pi; 0];

% Specify a time interval from 0 s to 10 s for finding the solution.
tInit  = 0;
tFinal = 10;

% Solve the ODE.
sols = ode45(f,[tInit tFinal],x0)

% Plot
figure;

yyaxis left;
plot(sols.x, sols.y(1,:), '-o');
ylabel('\theta (rad)');

yyaxis right;
plot(sols.x, sols.y(2,:), '-o');
ylabel('\theta_t (rad/s)');

grid on;
title('Nonlinear Motion using ode45 Solver');
xlabel('t (s)');