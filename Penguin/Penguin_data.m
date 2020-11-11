%% Pre-processor Stage 2
% Admin
clear;
% Import Raw Data
data = readmatrix('raw_data.txt');
%%
% Create array for mass meassurements from weigh-bridge
mass_data = data';
n_mass = 1:1:length(mass_data);

% Export to text file
writematrix(mass_data,'mass_measurements','FileType','text')
writematrix(n_mass,'mass_sample_size','FileType','text')

% Sampling rate of sensor
f = 50; % Hertz
dt = 1/f;

% Find gradient of the mass
gradient_mass_data = gradient(mass_data,dt);
n_gradient_mass = 1:1:length(gradient_mass_data);

% Export to etxt file
writematrix(gradient_mass_data,'gradient_mass_measurements','FileType','text')
writematrix(n_gradient_mass,'gradient_mass_sample_size','FileType','text')

%% Plotting
figure(1),clf
tiled = tiledlayout(2,1);

nexttile
plot(n_mass,mass_data,'-')
grid on
title('Mass Measurements from Weigh-bridge')
ylabel('Mass (kg)')
xlabel('Number of Samples')

nexttile
stairs(n_gradient_mass,gradient_mass_data,'Color','#D95319')
grid on
title('Gradient of Mass Measurements')
ylabel('Change in Mass (kg/s)')
xlabel('Number of Samples')