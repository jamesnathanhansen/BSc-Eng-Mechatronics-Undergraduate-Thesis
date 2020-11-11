%% Pre-proccessor Stage 3
% Admin
clear;

%% Preliminary Information 
% Import Data
n_mass = readmatrix('mass_sample_size.txt');
mass_data = readmatrix('mass_measurements.txt');
 
n_gradient_mass = readmatrix('gradient_mass_sample_size.txt');
gradient_mass_data = readmatrix('gradient_mass_measurements.txt');

%% Isolate Penguin Crossing Events
% thresholds
threshold_mass = 1.4;
gradient_threshold_mass = 15;

% Find start and end indices of each event within dataset 
reference = '';
run_index = {};

% Filter mass data
for i = 1:length(gradient_mass_data)
    if (mass_data(i)>threshold_mass & gradient_mass_data(i)>gradient_threshold_mass) & ~strcmp(reference,'start')
        reference = 'start';
        run_index = [run_index; {reference n_gradient_mass(i)}];
    elseif (mass_data(i)<threshold_mass & gradient_mass_data(i)<-gradient_threshold_mass) & ~strcmp(reference,'end')
        reference = 'end';
        run_index = [run_index; {reference n_gradient_mass(i)}];
    end
end

% Remove case where end is at the front of cell array
if strcmp(run_index(1,1),'end')
    run_index(1,:) = [];
end

% Remove case where start is at the end of cell array
if strcmp(run_index(end,1),'start')
    run_index(end,:) = [];
end

% Create cell array where each cell is one penguin crossing event
for i = 1:2:length(run_index)
    run = create_dataset(i,run_index,mass_data);
    run_set{i,1} = run;
end 

% Remove empty arrays
run_total = {};
for i = 1:2:length(run_set)
    run_total = [run_total; run_set{i}];
end 

% Append run time for each crossing to cell array
for i = 1:length(run_total)
    n_run = 1:length(run_total{i});
    run_total{i,2} = n_run;
end

%% Filter Penguin Crossing Events for KF measurement inputs
% percentage_cutoff Cutoff
percentage_cutoff = 0.05;

processed_run_total = {};
writecell(processed_run_total,'processed_run_total','FileType','text','Delimiter',',')

% Remove percentage_cutoff samples from start and end of each crossing event
for i = 1:length(run_total);
    temp = length(run_total{i,1});
    start = round(percentage_cutoff*temp);
    if start == 0
        start = 1;
    end
    finish = temp - round(percentage_cutoff*temp);
    processed_run_total{i,1} = run_total{i,1}(start:finish);
    processed_run_total{i,2} = run_total{i,2}(start:finish);
end

% Write to text file
for i = 1:length(processed_run_total)
    writecell(processed_run_total(i,1),'processed_run_total','FileType','text','WriteMode','append','Delimiter',',')
end
%% Write paramaters of threshold and percentage_cutoff cut-off to text file
parameters = [threshold_mass gradient_threshold_mass percentage_cutoff];
% writematrix(parameters,'changing_gradient_mass_threshold_parameters','FileType','text','WriteMode','append')

%% Measurement Noise Covariance, R
measurement_noise_mass = [];

% Isolate Data without Penguin Crossing Events (Sensor Noise)
for i = 1:length(mass_data)
    if mass_data(i)<0.3 && mass_data(i)>-0.4
        measurement_noise_mass = [measurement_noise_mass mass_data(i)];
    end
end 

measurement_noise_mass_n = 1:1:length(measurement_noise_mass);

% Find Measurment Noise Covariance
R_cov = cov(measurement_noise_mass);
% Export to text file
writematrix(R_cov,'R_cov','FileType','text')

%% Initial States x
x_int = {};
for i = 1:length(processed_run_total)
    x_int_mass = mean(processed_run_total{i,1});
    x_int{i,1} = x_int_mass;
end
% Export to text file
writecell(x_int,'x_int','FileType','text')

%% State Covariance, P
P_int = {};
for i = 1:length(processed_run_total)
%     P_int_mass = abs(x_int{i}-mean(processed_run_total{i,1}));
    P_int_mass = R_cov*10;
    P_int{i,1} = P_int_mass;
end
% Export to text file
writecell(P_int,'P_int','FileType','text')

%% Process Covariance Q
Q = {};
for i = 1:length(processed_run_total)
    Q_cov = cov(processed_run_total{i,1});
    spectral_density = 1000;
    Q{i,1} = Q_cov./spectral_density;
end
% Export to text file
writecell(Q,'Q_cov','FileType','text')

%% Plotting
figure(2),clf
tiled = tiledlayout('flow');
title(tiled, 'Penguin Crossing Events [cut-off = 5%] [thresholds: 1.4kg, 15kg/s]')
xlabel(tiled, 'Number of Samples')
ylabel(tiled, 'mass (kg)')
for i = 1:length(processed_run_total)
    nexttile
    hold on
    plot(run_total{i,2}, run_total{i,1})
    plot(processed_run_total{i,2}, processed_run_total{i,1},'LineWidth',1)
    title(i)
    grid on
    hold off
end

%% Function
% Function to append penguin crossing event to cell array
function [dataset] = create_dataset(position,data_index,data)
    dataset = data(data_index{position,2}:data_index{position+1,2});
end