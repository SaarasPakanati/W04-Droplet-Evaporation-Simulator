%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%           WD04: Droplet Evaporation Physical Setup Simulation           %
%                                                                         %
%                                                                         %
% This code simulates a small subset of the planned experiments on the    %
%   the droplet evaporation experiments conducted by the Senior Design    %
%   team.                                                                 %
%                                                                         %
% Developed By WD04 Droplet Evaporation Senior Design Team 2025-26,       %
% The University of Cincinnati, February 16th 2025                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Initialization
clc; clear; close all;

%% Model Paramaters

% Physical Parameters [CONFIRM THESE PARAMETERS!]
T_s = 80;                                       % [K] Substrate Temperature
V_uL = 300;                                     % [uL] Initial Droplet Volume
t = 100;%500;                                        % [s] Length of Simulation
dt = 1;                                         % [s] Temporal Discretization

% Physical Parameters [YOU CAN CHANGE, BUT DO YOU REALLY WANT TO?]
k      = 0.6515;                                % [W/m-K] Thermal Conductivity
rho_l  = 983.2;                                 % [kg/m^3] Liquid Density
sigma  = 0.06624;                               % [N/m] Surface Tension
hfg    = 2358e3;                                % [J/kg] Latent Heat of Vaporization

% Physical Parameters [DO NOT CHANGE!]
h_conv = 1e+2;                                  % [W/m^2-K] Convection Coefficient
alpha  = 1.1e-9;                                % [-] Accommodation Coefficient
beta   = 0.975;                                 % [-] Correction Factor
theta_deg = 10;                                 % [deg] Contact Angle

% PID CONTROL PARAMETERS [DO NOT CHANGE]
use_PID = true;                                 % Toggle Control ON/OFF
Kp = 1e5;                                       % Proportional Gain
Ki = 5e3;                                       % Integral Gain
Kd = 0;                                         % Derivative Gain
pid_integral = 0;                               % Accumulator initialization
pid_prev_error = 0;                             % Derivative initialization


% Array Initialization
massFlux = nan(size(1:dt:t));
xValues_last = nan(size(1:dt:t));
xValues_star_last = nan(size(1:dt:t));
hValues = nan(size(1:dt:t));
hValues_star = nan(size(1:dt:t));
Volume = nan(size(1:dt:t));
hVideo = nan(3e5, length(1:dt:t));
xVideo = nan(3e5, length(1:dt:t));

% Calculate Initial Radius (For PID control)
[~, ~, x_init, ~] = volumeToGeometry(V_uL, theta_deg);
targetRadius = max(x_init); 
fprintf('Target Contact Radius: %.6f mm\n', targetRadius*1e3);


for i = 1:length(1:dt:t)

    % Numerical Model
    [m_total, Ti_avg, x, h, converged] = evaporationModel( V_uL, T_s, h_conv, alpha, beta, theta_deg, k, rho_l, sigma, hfg);
    
    V_loss = (m_total * dt / rho_l) * 1e9;                                      % [uL] Volume lost
    
    V_temp = V_uL + V_loss;                                                     % [uL] Current Droplet Volume
    
    % Geometry of Droplet after lost liquid due to evaporation
    [~, ~, x_temp, ~] = volumeToGeometry(V_temp, theta_deg);
    currentRadius = max(x_temp);
    
    % PID Controller
    V_add = 0;
    if use_PID
        [V_add, pid_integral, pid_prev_error] = pidVolumeControl(targetRadius, currentRadius, dt, Kp, Ki, Kd, pid_integral, pid_prev_error);
    end
    
    % Final Volume after PID addition.
    V_uL = V_temp + V_add;
    
    % Calculate height profile after evaporated mass and PID control
    [h_apex, theta_c, x_star, h_star] = volumeToGeometry(V_uL, theta_deg);

    % Adding values
    massFlux(i) = m_total;                      % [kg/s] Mass Evaporation Rate
    xValues_last(i) = max(x);                   % [m] Location of Contact Line
    xValues_star_last(i) = max(x_star);         % [m] NEW Location of Contact Line
    hValues(i) = max(h);                        % [m] Height Profile
    hValues_star(i) = max(h_star);              % [m] NEW Height Profile
    Volume(i) = V_uL;                           % [uL] Volume of Droplet
    hVideo(1:length(h), i) = h';
    xVideo(1:length(x), i) = x';

end

% Plots
close all;
figure(Position=[100 100 1200 700]);

% Total Evaporation Rate over Time
subplot(2, 3, 1);
plot(1:dt:t, -massFlux*1e6, 'b-', 'LineWidth', 2);
xlabel('Time [s]');
ylabel('Rate [mg/s]');
title('Evaporation Rate');
grid on;

% Contact Line Position over Time
subplot(2, 3, 2);
plot(1:dt:t, xValues_last * 1e3, 'r-', 'LineWidth', 2);
hold on;

% Safety check: If target_contact_radius is missing for some reason, use the first value
if ~exist('target_contact_radius', 'var')
    target_contact_radius = xValues_last(1);
end

yline(target_contact_radius*1e3, 'k--', 'Target');
xlabel('Time [s]');
ylabel('Position [mm]');
title('Contact Line Position');
legend('Actual', 'Target', 'Location','southwest');
grid on;
if ~all(isnan(xValues_last))
    ylim([min(xValues_last*1e3)*0.99, max(xValues_last*1e3)*1.01]);
end

% Apex Height over Time
subplot(2, 3, 3);
plot(1:dt:t, hValues * 1e3, 'g-', 'LineWidth', 2);
xlabel('Time [s]');
ylabel('Height [mm]');
title('Apex Height');
grid on;

% Volume over Time
subplot(2, 3, 4);
plot(1:dt:t, Volume, 'k-', 'LineWidth', 2);
xlabel('Time [s]');
ylabel('Volume [uL]');
title('Droplet Volume');
grid on;

% Error
subplot(2, 3, 5);
error_mm = (target_contact_radius - xValues_last) .* 1e3;
plot(1:dt:t, error_mm, 'c-', 'LineWidth', 2);
xlabel('Time [s]');
ylabel('Error [mm]');
title('Control Error');
grid on;

sgtitle('Droplet Evaporation with PID Level Control');

%%
figure(Position=[488 241.8000 560 158.4000]);
tArray = 1:dt:t;

v = VideoWriter("DropletSizePID_1.avi");
open(v)

for i = 1 : 5 : length(hVideo(1, :))
    plot(rmmissing(xVideo(:, i))*1e3, rmmissing(hVideo(:, i))*1e3, 'c-', 'LineWidth', 5);
    xlabel('x [mm]');
    ylabel('h [mm]');
    title(sprintf("Height Profile @ t = %.2f", tArray(i)));
    % grid on;
    ylim([0 1]);
    xlim([-10 10]);
    % ylim([-0.001 0.02]);
    % xlim([11.3 11.35]);
    drawnow();
    frame = getframe(gcf);
    writeVideo(v,frame);
end

close(v);
    