% evaporationModel Function
%       This function will calculate the total evaporation rate and average
%       interface temperature.
%
% [m_total, Ti_avg, h, converged] = evaporationModel( V_uL, T_s, h_conv, alpha, beta, theta_deg, k, rho_l, sigma, hfg)
%
% Input(s):
%   m_total ------------------------ [float] Integrated Evaporation Rate
%   Ti_avg ------------------------- [float] Averaged Interface Temperature
%   x ------------------------------ [1-D] Spatial Mapping
%   h ------------------------------ [1-D] Height Profile
%   converged ---------------------- [boolean] Convergence Check
%
% Output(s):
%   V_uL --------------------------- [float] Volume of the droplet
%   T_s ---------------------------- [float] Substrate Temperature
%   h_conv ------------------------- [float] Convection Coefficient
%   alpha -------------------------- [float] Accomodation Coefficient
%   beta --------------------------- [float] Beta Factor
%   theta_deg ---------------------- [float] 3 Phase Contact Angle of the Droplet
%   rho_l -------------------------- [float] Liquid Density
%   sigma -------------------------- [float] Surface Tension
%   hfg ---------------------------- [float] Latent heat of Evaporation
%
% Developed By: Saaras Pakanati @ The University of Cincinnati, February 16 2026

function [m_total, Ti_avg, x, h, converged] = evaporationModel( V_uL, T_s, h_conv, alpha, beta, theta_deg, k, rho_l, sigma, hfg)

    % Universal constants
    R_g = 8.314462618153240;                                        % [J/K-mol] Ideal Gas Constant
    M   = 18.015e-3;                                                % [kg/mol] Molecular Mass of H2O
    T_inf = 25;                                                     % [C] Room Temperature
    T_inf_K = T_inf + 273.15;                                       % [K] Room Temperature
    T_s_K = T_s + 273.15;                                           % [K] Surface Temperature
    
    % Geometry from volume
    [h_apex, theta_c] = volumeToGeometry(V_uL, theta_deg);
    R = h_apex / (1 - cos(theta_c));                                % [m] Radius of Curvature
    y_c = -R + h_apex;                                              % [m] Circle Center Correction
    a = R * sin(theta_c);                                           % [m] Contact Line radius
    
    % Discretization
    dx = 1e-7;                                                      % [m] Spatial Discretization
    x = -a : dx : a;                                                % Spatial Domain Generation
    h = y_c + sqrt(R^2 - x.^2);                                     % [m] Droplet Height Profile
    h(h < 0) = 0;
    
    % Saturation pressure (Antoine)
    A = 8.07131; B = 1730.63; C = 233.426;
    pv_sat = 10^(A - B/(C + T_s)) * 133.322;                        % [Pa] Saturation Pressure
    pv_inf = 10^(A - B/(C + T_inf)) * 133.322;                      % [Pa] Ambient Pressure
    pv = pv_inf;                                                    % [Pa] Ambient Vapor Pressure
    
    % Vapor properties
    Tv = T_s_K;                                                     % [K] Vapor Temperature
    rho_v = (pv * M) / (R_g * Tv);                                  % [kg/m^3] Vapor Density
    
    % Initial guess for Interface Iemperature (No Evaporation)
    Ti = T_s_K - (h_conv .* h .* (T_s_K - T_inf_K)) ./ (k + (h_conv .* h));
    Ti_old = Ti;
    
    % Constant Curvature (Spherical Cap Assumption)
    kappa = 1 / R;
    
    % Itterative solver parameters
    max_iter = 1E+6;                                                % Maximum Allowable Itterations
    tol = 1e-6;                                                     % Convergence Criterion
    relax = 0.5;                                                    % Relaxation Factor
    converged = false;                                              % Converged Boolean
    
    for iter = 1:max_iter
        % Evaporation Model (Bellur)
        W1 = pv_sat ./ pv;
        W2 = (1 - (Tv ./ Ti)) .* ((pv_sat .* hfg) ./ pv);
        W3 = (Tv ./ Ti) .* (rho_v ./ rho_l) .* (sigma .* kappa ./ pv);
        W = W1 + W2 + W3;
        
        bellur_1 = (2 .* alpha) ./ (2 - alpha);
        bellur_2 = sqrt(M ./ (2 .* pi .* R_g .* Tv));
        bellur_3 = pv .* beta .* W .* sqrt(Tv ./ Ti) - 1;
        m = bellur_1 .* bellur_2 .* bellur_3;                       % [kg/m^2-s] Local Mass Flux
        
        % Energy balance at interface: Conduction = Convection + Evaporation
        Ti_new = (k .* T_s_K ./ h + h_conv .* T_inf_K - m .* hfg) ./ (k ./ h + h_conv);
        
        Ti = relax .* Ti_new + (1 - relax) .* Ti_old;               % Relaxation for stability
        
        % Convergence Check
        max_change = max(abs(Ti - Ti_old));
        if max_change < tol
            converged = true;
            break;
        end
        
        Ti_old = Ti;                                                % Update Temperature
    end
    
    % Integrate total evaporation rate over surface (axisymmetric)
    m_total = 0;  % [kg/s]
    for idx = 2 : length(x)
        dh_dx = (h(idx) - h(idx-1)) / dx;                           % Slope
        dS_factor = sqrt(1 + dh_dx^2);                              % Area Correction
        r = abs(x(idx));                                            % Radial Coordinate
        m_total = m_total + 2 * pi * r * m(idx) * dS_factor * dx;
    end
    
    Ti_avg = mean(Ti);                                              % [K] Average Interface Temperature
end