clear;

performance_and_parameters = readmatrix('reduced_runs_performance.txt');

figure
title('Std. Dev. & Mean vs Minimum Number of Samples [cut-off = 5%, thresholds: 1.4kg, 15kg/s]')
xlabel('Minimum Number of Samples')
grid on

yyaxis left;
plot(performance_and_parameters(:,1), performance_and_parameters(:,2),'-o')
ylabel('Std. Dev \sigma')

yyaxis right;
plot(performance_and_parameters(:,1), performance_and_parameters(:,3),'-o')
ylabel('Average \mu')