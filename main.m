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

% Model Paramaters

% Physical Parameters [CONFIRM THESE PARAMETERS!]
T_s = 90;                                       % [K] Substrate Temperature
V_uL = 200;                                     % [uL] Initial Droplet Volume
t = 300;                                        % [s] Length of Simulation
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

% Array Initialization
massFlux = nan(size(1:dt:t));
xValues_last = nan(size(1:dt:t));
xValues_star_last = nan(size(1:dt:t));
hValues = nan(size(1:dt:t));
hValues_star = nan(size(1:dt:t));
Volume = nan(size(1:dt:t));

for i = 1:length(1:dt:t)

    % Numerical Model
    [m_total, Ti_avg, x, h, converged] = evaporationModel( V_uL, T_s, h_conv, alpha, beta, theta_deg, k, rho_l, sigma, hfg);
    
    V_uL_star = V_uL + (m_total * dt / rho_l) * 1e9;            % [uL] Volume remaining after evaporation
    
    % Calculate height profile after evaporated mass.
    [h_apex, theta_c, x_star, h_star] = volumeToGeometry(V_uL_star, theta_deg);

    % Adding values
    massFlux(i) = m_total;                      % [kg/s] Mass Evaporation Rate
    xValues_last(i) = max(x);                   % [m] Location of Contact Line
    xValues_star_last(i) = max(x_star);         % [m] NEW Location of Contact Line
    hValues(i) = max(h);                        % [m] Height Profile
    hValues_star(i) = max(h_star);              % [m] NEW Height Profile
    Volume(i) = V_uL;                           % [uL] Volume of Droplet
    
    % Update the volume of the droplet.
    V_uL = V_uL_star;

end

%% Plots
close all;
figure(Position=[488 87 871 574]);

% Total Evaporation Rate over Time
subplot(2, 3, 1);
plot(1:dt:t, -massFlux*1e6, 'b-', 'LineWidth', 2);
xlabel('Time [s]');
ylabel('Mass Evaporation Rate [mg/s]');
title('Total Evaporation Rate');
grid on;

% Contact Line Position over Time
subplot(2, 3, 2);
plot(1:dt:t, xValues_last * 1e3, 'r-', 'LineWidth', 2);
xlabel('Time [s]');
ylabel('Contact Line Position [mm]');
title('Contact Line Position');
grid on;

% Apex Height over Time
subplot(2, 3, 3);
plot(1:dt:t, hValues * 1e3, 'g-', 'LineWidth', 2);
xlabel('Time [s]');
ylabel('Apex Height [mm]');
title('Apex Height');
grid on;

% Volume over Time
subplot(2, 3, 4);
plot(1:dt:t, Volume, 'k-', 'LineWidth', 2);
xlabel('Time [s]');
ylabel('Volume [uL]');
title('Droplet Volume');
grid on;

% Change in Contact Line position
subplot(2, 3, 5);
delta_x = (xValues_last - xValues_last(1)) * 1e3;
plot(1:dt:t, delta_x, 'm-', 'LineWidth', 2);
xlabel('Time [s]');
ylabel('delta x [mm]');
title('Change in Contact Line Position');
grid on;

% Change in Apex Height
subplot(2, 3, 6);
delta_h = (hValues - hValues(1)) * 1e3;
plot(1:dt:t, delta_h, 'c-', 'LineWidth', 2);
xlabel('Time [s]');
ylabel('delta h mm]');
title('Change in Apex Height');
grid on;

sgtitle('Droplet Evaporation Dynamics');