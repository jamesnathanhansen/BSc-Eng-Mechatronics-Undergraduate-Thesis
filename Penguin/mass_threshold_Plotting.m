clear;

parameters = readmatrix('changing_mass_threshold_parameters.txt');
performance = readmatrix('changing_mass_threshold_performance.txt');

figure
title('Std. Dev. & Mean vs Mass Threshold [Cut-off = 5%, Gradient of Mass Threshold = 10 kg/s]')
xlabel('Mass Threshold (kg)')
grid on

yyaxis left;
plot(parameters(:,1), performance(:,1),'-o')
ylabel('Std. Dev \sigma')

yyaxis right;
plot(parameters(:,1), performance(:,2),'-o')
ylabel('Average \mu')