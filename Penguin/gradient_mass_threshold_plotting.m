clear;

parameters = readmatrix('changing_gradient_mass_threshold_parameters.txt');
performance = readmatrix('changing_gradient_mass_threshold_performance.txt');

figure
title('Std. Dev. & Mean vs Gradient of Mass Threshold [Cut-off = 5%, Mass Threshold = 1.4kg]')
xlabel('Graident of Mass Threshold (kg/s)')
grid on

yyaxis left;
plot(parameters(:,2), performance(:,1),'-o')
ylabel('Std. Dev \sigma')

yyaxis right;
plot(parameters(:,2), performance(:,2),'-o')
ylabel('Average \mu')