% volumeToGeometry Function
%       This function will calculate the apex height and contact angle to
%       compute the height profile of the droplet.
%
% [h_apex, theta_c, x_star, h_star] = volumeToGeometry(V_uL, theta_deg)
%
% Input(s):
%   V_uL --------------------------- [float] Volume of Droplet
%   theta_deg ---------------------- [float] Contact Angle
%   x_star ------------------------- [1-D] Spatial Mapping
%   h_star ------------------------- [1-D] Height Profile
%
% Output(s):
%   h_apex ------------------------- [float] Apex Height of Droplet
%   theta_c ------------------------ [float] 3 Phase Contact Angle of the Droplet
%
% Developed By: Saaras Pakanati @ The University of Cincinnati, February 16 2026

function [h_apex, theta_c, x_star, h_star] = volumeToGeometry(V_uL, theta_deg)
    
    theta_c = theta_deg * pi / 180;                                 % [rad] Contact Angle
    V = V_uL * 1e-9;                                                % [m^3] Volume
    
    term = (2 + cos(theta_c)) / (1 - cos(theta_c));   
    h_apex = (3 * V / (pi * term))^(1/3);                           % [m] Apex Height


    R = h_apex / (1 - cos(theta_c));                                % [m] Radius of Curvature
    y_c = -R + h_apex;                                              % [m] Circle Center Correction
    a = R * sin(theta_c);                                           % [m] Contact Line radius
    
    % Discretization
    dx = 1e-7;                                                      % [m] Spatial Discretization
    x = -a : dx : a;                                                % Spatial Domain Generation
    h = y_c + sqrt(R^2 - x.^2);                                     % [m] Droplet Height Profile
    h(h < 0) = 0;
    x_star = x;
    h_star = h;
end