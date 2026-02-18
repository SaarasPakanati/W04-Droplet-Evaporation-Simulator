% pidVolumeControl Function
%       This function calculates the necessary volume addition to maintain
%       a target contact line position using a PID controller.
%
% [V_add, integral_error, prev_error] = pidVolumeControl(target_pos, current_pos, dt, Kp, Ki, Kd, integral_error, prev_error)
%
% Input(s):
%   targetPos ---------------------- [float] Target Contact Line Position
%   currentPos --------------------- [float] Current Contact Line Position
%   dt ----------------------------- [float] Time Step
%   Kp ----------------------------- [float] P Gain
%   Ki ----------------------------- [float] I Gain
%   Kd ----------------------------- [float] D Gain
%   integralError ------------------ [float] Accumulated Error
%   prevError ---------------------- [float] Previous Temporal Step Error
%
% Output(s):
%   V_add -------------------------- [float] Volume to add
%   integralError ------------------ [float] Updated Accumulated Error
%   prevError ---------------------- [float] Updated Previous Error

function [V_add, integralError, prevError] = pidVolumeControl(targetPos, currentPos, dt, Kp, Ki, Kd, integralError, prevError)
    
    error = targetPos - currentPos;                             % Check if droplet is shrinking and error
    
    % PID Terms Comutation
    P = Kp * error;                                             % Proportional Term
    integralError = integralError + (error * dt);
    I = Ki * integralError;                                     % Integral Term
    derivative = (error - prevError) / dt;
    D = Kd * derivative;                                        % Derivative Term
    
    % Computing output volume
    u = P + I + D;
    
    % ONLY ADD liquid, never take out. (APhysical for the setup...)
    if u < 0
        V_add = 0;
    else
        V_add = u;
    end
    
    prevError = error;                                          % Store Error

end