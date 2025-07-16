%========================================================================%
% This function calculates the forces and the torque caused to the object by
% the ground, given:
%     - mass
%     - mu (coefficient of friction)
%     - contact_area (are of the object in touch with the ground)
%     - ground_friction_parameter (simulation related parameter, tunable)
%     - dp (slight change of the position and orientation)
%========================================================================%

function [friction, number] = calculate_friction_with_ground(L, dp, ground_friction_parameter)

    if sqrt(dp' * inv(L) * dp) > ground_friction_parameter
        number = 1;
        friction = inv(L) * dp / sqrt(dp' * inv(L) * dp);
    else
        number = 0;
        friction = inv(L) * dp / ground_friction_parameter;
    end

end