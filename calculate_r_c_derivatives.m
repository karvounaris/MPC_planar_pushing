function r_c_partial_derivative_phi = calculate_r_c_derivatives(phi, radius, len, width, object_shape)
%=======================================================================%
% Calculates the partial derivatives (dx_c/dphi, dy_c/dphi) of the
% contact point (x_c, y_c) in the object body frame, for two shapes:
%   1) "rectangular_capsule_prism"
%   2) "rectangular_prism"
%
% Inputs:
%   phi          : angle [0..2*pi)
%   radius       : radius used by the capsule shape
%   len, width   : rectangle dimensions (for the rectangular prism)
%   object_shape : string specifying shape type
%
% Output:
%   r_c_partial_derivative_phi = [ dx_c/dphi; dy_c/dphi ]
%=======================================================================%

    % Normalize phi to [0, 2*pi)
    phi = mod(phi, 2*pi);
    
    if object_shape == "rectangular_capsule_prism"

        if phi >= 0 && phi < pi/4
            x_c_partial_derivative_phi = -2 * radius * sin(2*phi);
            y_c_partial_derivative_phi =  2 * radius * cos(2*phi);

        elseif phi >= pi/4 && phi < 3*pi/4
            x_c_partial_derivative_phi = -radius * csc(phi)^2; 
            y_c_partial_derivative_phi =  0;

        elseif phi >= 3*pi/4 && phi < 5*pi/4
            x_c_partial_derivative_phi =  2 * radius * sin(2*phi);
            y_c_partial_derivative_phi = -2 * radius * cos(2*phi);

        elseif phi >= 5*pi/4 && phi < 7*pi/4
            x_c_partial_derivative_phi =  radius * csc(phi)^2;
            y_c_partial_derivative_phi =  0;

        elseif phi >= 7*pi/4 && phi < 2*pi
            x_c_partial_derivative_phi = -2 * radius * sin(2*phi);
            y_c_partial_derivative_phi =  2 * radius * cos(2*phi);

        else
            x_c_partial_derivative_phi = 0;
            y_c_partial_derivative_phi = 0;
        end

    elseif object_shape == "rectangular_prism"

        alpha = atan(width/len);

        if ((phi >= 0) && (phi < alpha)) || ((phi >= 2*pi - alpha) && (phi < 2*pi))
            % Right edge:
            %   x_c = +len/2        => dx/dphi = 0
            %   y_c = (len/2)*tan(phi)
            %   => dy/dphi = (len/2)*sec^2(phi)
            x_c_partial_derivative_phi = 0;
            y_c_partial_derivative_phi = (len/2)*sec(phi)^2;

        elseif (phi >= alpha) && (phi < pi - alpha)
            % Top edge:
            %   x_c = (width/2)*cot(phi) = (width/2)*cos(phi)/sin(phi)
            %       => dx/dphi = -(width/2)*csc^2(phi)
            %   y_c = +width/2
            %       => dy/dphi = 0
            x_c_partial_derivative_phi = -(width/2)*csc(phi)^2;
            y_c_partial_derivative_phi = 0;

        elseif (phi >= pi - alpha) && (phi < pi + alpha)
            % Left edge:
            %   x_c = -len/2        => dx/dphi = 0
            %   y_c = -(len/2)*tan(phi)
            %       => dy/dphi = -(len/2)*sec^2(phi)
            x_c_partial_derivative_phi = 0;
            y_c_partial_derivative_phi = -(len/2)*sec(phi)^2;

        elseif (phi >= pi + alpha) && (phi < 2*pi - alpha)
            % Bottom edge:
            %   x_c = - (width/2)*cot(phi)
            %       => dx/dphi = - (width/2)* d/dphi[cot(phi)]
            %                    =  (width/2)*csc^2(phi)
            %   y_c = -width/2
            %       => dy/dphi = 0
            x_c_partial_derivative_phi = +(width/2)*csc(phi)^2;
            y_c_partial_derivative_phi = 0;

        else
            % Possibly exactly at a corner angle or invalid phi
            x_c_partial_derivative_phi = 0;
            y_c_partial_derivative_phi = 0;
        end

    else
        % In case object_shape is none of the above
        x_c_partial_derivative_phi = 0;
        y_c_partial_derivative_phi = 0;
    end

    % Combine as a 2x1 vector
    r_c_partial_derivative_phi = [x_c_partial_derivative_phi; y_c_partial_derivative_phi];
end
