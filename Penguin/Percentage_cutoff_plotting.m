clear;

performance = readmatrix('changing_percentage_cutoff_performance.txt');
parameters = readmatrix('changing_percentage_cutoff_parameters.txt');

figure
title('Std. Dev. & Mean vs Percentage Cut-off [thresholds: 1kg, 10kg/s]')
xlabel('% Cut-off')
grid on

yyaxis left;
plot(parameters(:,3).*100, performance(:,1),'-o')
ylabel('Std. Dev \sigma')

yyaxis right;
plot(parameters(:,3).*100, performance(:,2),'-o')
ylabel('Average \mu')