function [x_c, y_c, n_c, t_c] = calculate_r_c_casadi(phi, len, radius, wid, object_shape)
    import casadi.*

    % Wrap phi into [0, 2*pi)
    phi_mod = fmod(phi, 2*pi);
    phi_mod = if_else(phi_mod < 0, phi_mod + 2*pi, phi_mod);

    % --------------------------------------------------------
    if object_shape == "rectangular_capsule_prism"
        % === Define conditions ===
        cond1 = (phi_mod >= 0) & (phi_mod < pi/4);
        cond2 = (phi_mod >= pi/4) & (phi_mod < 3*pi/4);
        cond3 = (phi_mod >= 3*pi/4) & (phi_mod < 5*pi/4);
        cond4 = (phi_mod >= 5*pi/4) & (phi_mod < 7*pi/4);
        % fallback is [7pi/4, 2pi)

        % === Define expressions ===
        x_c1 = len/2 + radius * cos(2*phi_mod);
        y_c1 = radius * sin(2*phi_mod);
        n_c1 = [-cos(2*phi_mod); -sin(2*phi_mod)];
        t_c1 = [sin(2*phi_mod); -cos(2*phi_mod)];

        x_c2 = radius * (cos(phi_mod)/sin(phi_mod));
        y_c2 = radius;
        n_c2 = [0; -1];
        t_c2 = [1; 0];

        x_c3 = -(len/2 + radius * cos(2*phi_mod));
        y_c3 = -radius * sin(2*phi_mod);
        n_c3 = [cos(2*phi_mod); sin(2*phi_mod)];
        t_c3 = [-sin(2*phi_mod); cos(2*phi_mod)];

        x_c4 = -radius * (cos(phi_mod)/sin(phi_mod));
        y_c4 = -radius;
        n_c4 = [0; 1];
        t_c4 = [-1; 0];

        x_c5 = x_c1;
        y_c5 = y_c1;
        n_c5 = n_c1;
        t_c5 = t_c1;

        % === Combine with if_else ===
        x_c = if_else(cond1, x_c1, ...
               if_else(cond2, x_c2, ...
               if_else(cond3, x_c3, ...
               if_else(cond4, x_c4, x_c5))));

        y_c = if_else(cond1, y_c1, ...
               if_else(cond2, y_c2, ...
               if_else(cond3, y_c3, ...
               if_else(cond4, y_c4, y_c5))));

        n_c = if_else(cond1, n_c1, ...
               if_else(cond2, n_c2, ...
               if_else(cond3, n_c3, ...
               if_else(cond4, n_c4, n_c5))));

        t_c = if_else(cond1, t_c1, ...
               if_else(cond2, t_c2, ...
               if_else(cond3, t_c3, ...
               if_else(cond4, t_c4, t_c5))));

    % --------------------------------------------------------
    elseif object_shape == "rectangular_prism"
        % === Rectangular prism ===
        alpha = atan(wid/len);

        cond1 = ((phi_mod >= 0) & (phi_mod < alpha)) | ((phi_mod >= 2*pi - alpha) & (phi_mod < 2*pi));
        cond2 = (phi_mod >= alpha) & (phi_mod < pi - alpha);
        cond3 = (phi_mod >= pi - alpha) & (phi_mod < pi + alpha);
        cond4 = (phi_mod >= pi + alpha) & (phi_mod < 2*pi - alpha);

        % Right edge
        x_c1 = len/2;
        y_c1 = (len/2)*tan(phi_mod);
        n_c1 = [-1; 0];
        t_c1 = [0; -1];

        % Top edge
        x_c2 = wid/2 * (cos(phi_mod)/sin(phi_mod));
        y_c2 = wid/2;
        n_c2 = [0; -1];
        t_c2 = [1; 0];

        % Left edge
        x_c3 = -len/2;
        y_c3 = -len/2 * tan(phi_mod);
        n_c3 = [1; 0];
        t_c3 = [0; 1];

        % Bottom edge
        x_c4 = -(wid/2)*cos(phi_mod)/sin(phi_mod);
        y_c4 = -wid/2;
        n_c4 = [0; 1];
        t_c4 = [-1; 0];

        % fallback: use right edge again
        x_c5 = x_c1;
        y_c5 = y_c1;
        n_c5 = n_c1;
        t_c5 = t_c1;

        % === Combine ===
        x_c = if_else(cond1, x_c1, ...
               if_else(cond2, x_c2, ...
               if_else(cond3, x_c3, ...
               if_else(cond4, x_c4, x_c5))));

        y_c = if_else(cond1, y_c1, ...
               if_else(cond2, y_c2, ...
               if_else(cond3, y_c3, ...
               if_else(cond4, y_c4, y_c5))));

        n_c = if_else(cond1, n_c1, ...
               if_else(cond2, n_c2, ...
               if_else(cond3, n_c3, ...
               if_else(cond4, n_c4, n_c5))));

        t_c = if_else(cond1, t_c1, ...
               if_else(cond2, t_c2, ...
               if_else(cond3, t_c3, ...
               if_else(cond4, t_c4, t_c5))));

    else
        error('Unknown object_shape for calculate_r_c_casadi');
    end
end
