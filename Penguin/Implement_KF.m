%% Extended Kalman Filter Implementation
%% Admin
clear;
%% Import Paramters
Q_cov = readcell('Q_cov.txt');
R_cov = readmatrix('R_cov.txt');
P_int = readcell('P_int.txt');
x_int_mass = readcell('x_int.txt');
measurements = readcell('processed_run_total.txt','ConsecutiveDelimitersRule','join');

z_kg = cell(size(measurements,1),1);
for i = 1:size(measurements)
    temp = [];
    for j = 1:length(measurements)
        if ismissing(measurements{i,j}) == 1
            continue
        else
            temp = [temp measurements{i,j}];
        end 
    end
    z_kg{i,1} = temp;
end
%% Implement EKF
final_masses = zeros(length(z_kg),1);
z_out = cell(size(z_kg,1),1);
xp_out = cell(size(z_kg,1),1);
xs_out = cell(size(z_kg,1),1);
P = cell(size(P_int,1),1);
t = cell(size(z_kg,1),1);

for i = 1:length(z_kg)
    [ts,xs,xp,z,Pvals] = Kalman_filter_mass(Q_cov{i},R_cov,P_int{i},x_int_mass{i},z_kg{i});
    z_out{i,1} = z;
    xp_out{i,1} = xp;
    xs_out{i,1} = xs;
    t{i,1} = ts;
    P{i,1} = Pvals;
    
    final_masses(i,1) = xs(end);
end 

%% Calculate Performance Parameters
std_dev = std(final_masses);
ave = mean(final_masses);
performance = [std_dev,ave]

%% Post-process EKF Output
xs_processed = {};
z_processed = {};
tp = {};
final_masses_processed = [];
P_processed = {};

% Only include penguin crossing events that have 40 samples or more (0.8
% seconds or longer)
min_num_samples = 40;
for i = 1:length(xs_out)
    if length(xs_out{i,1}) >= min_num_samples
        xs_processed{end+1} = xs_out{i,1};
        z_processed{end+1} = z_out{i,1};
        tp{end+1} = t{i,1};
        P_processed{end+1} = P{i,1};
        
        final_masses_processed = [final_masses_processed; xs_out{i,1}(end)];
    end 
end 

%% Calculate Processed Performance Parameters
std_dev = std(final_masses_processed);
ave = mean(final_masses_processed);
processed_performance = [std_dev,ave];

LP = 25;
final_masses_processed_nooutliers = rmoutliers(final_masses_processed,'percentiles', [LP 100]);
std_dev = std(final_masses_processed_nooutliers);
ave = mean(final_masses_processed_nooutliers);
processed_performance_nooutliers = [LP,std_dev,ave,length(final_masses_processed_nooutliers)];

% writematrix(processed_performance_nooutliers,'remove_outliers_performance','FileType','text','WriteMode','append')
std_dev
ave
%% Plot
% Plot xp (Predicted Output of States w/o Updates)
figure(3),clf
tiled = tiledlayout('flow');
title(tiled, 'Predicted Output of States without EKF Updates vs Measured States')
xlabel(tiled, 't (s)')
ylabel(tiled, 'Mass (kg)')
for i = 1:length(xp_out)
    nexttile
    hold on
    plot(t{i,1}, xp_out{i,1})
    plot(t{i,1}(2:end), z_out{i,1})
    title(i)
    grid on
    hold off
end

% Plot EKF Estimated States
figure(4),clf
tiled = tiledlayout('flow');
title(tiled, 'EKF Estimated States vs Measured States [Q = Q/1000]')
xlabel(tiled, 't (s)')
ylabel(tiled, 'Mass (kg)')
for i = 1:length(xs_out)
    nexttile
    hold on
    plot(t{i,1}, xs_out{i,1})
    plot(t{i,1}(2:end), z_out{i,1})
    title(i)
    grid on
    hold off
end


% Plot P_cov (Error Covariances)
figure(5),clf
tiled = tiledlayout('flow');
title(tiled, 'State Error Covaraince')
xlabel(tiled, 't (s)')
ylabel(tiled, 'Variance')
for i = 1:length(P)
    nexttile
    plot(t{i,1}, P{i,1})
    title(i)
    grid on
end

% Plot Processed EKF Estimated States
figure(6),clf
tiled = tiledlayout('flow');
title(tiled, 'Processed EKF Estimated States vs Measured States [Q = Q/1000,  min. # of samples = 40]')
xlabel(tiled, 't (s)')
ylabel(tiled, 'Mass (kg)')
for i = 1:length(xs_processed)
    nexttile
    hold on
    plot(tp{1,i}, xs_processed{1,i})
    plot(tp{1,i}(2:end), z_processed{1,i})
    title(i)
    grid on
    hold off
end

% Plot Processed P_cov (Error Covariances)
figure(7),clf
tiled = tiledlayout('flow');
title(tiled, 'Processed State Error Covaraince')
xlabel(tiled, 't (s)')
ylabel(tiled, 'Variance')
for i = 1:length(P_processed)
    nexttile
    plot(tp{1,i}, P_processed{1,i})
    title(i)
    grid on
end