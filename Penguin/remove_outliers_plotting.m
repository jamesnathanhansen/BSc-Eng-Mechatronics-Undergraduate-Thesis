clear;

performance_and_parameters = readmatrix('remove_outliers_performance_and_param.txt');

figure
tiled = tiledlayout(2,1);
xlabel(tiled, 'Percentile (%)')
nexttile
title('Std. Dev. & Mean vs Lower Percentile Outlier Removal [cut-off = 5%, thresholds: 1.4kg 15kg/s, min # of samples = 40]')
grid on

yyaxis left;
plot(performance_and_parameters(:,1), performance_and_parameters(:,2),'-o')
ylabel('Std. Dev \sigma')

yyaxis right;
plot(performance_and_parameters(:,1), performance_and_parameters(:,3),'-o')
ylabel('Average \mu')

nexttile
plot(performance_and_parameters(:,1), performance_and_parameters(:,4),'-o')
ylabel('Number of Penguin Crossing Events')
grid on

