function [A_fun, B_fun] = casadi_linearize_system(L, radius, len, wid, object_shape)

    import casadi.*

    % === 1) Symbolic variables ===
    x_vec = MX.sym('x', 4);   % [x; y; theta; phi]
    u_vec = MX.sym('u', 3);   % [fn; ft; phi_dot]

    theta = x_vec(3);
    phi   = x_vec(4);

    fn = u_vec(1);
    ft = u_vec(2);
    phi_dot = u_vec(3);

    l1 = L(1,1);
    l2 = L(2,2);
    l3 = L(3,3);

    alpha = atan(wid/len);

    % === 2) Capsule piecewise cases ===
    dx_capsule1 = [
        (l1*(-cos(2*phi))*cos(theta) - l2*(-sin(2*phi))*sin(theta))*fn + ...
        (l1*(sin(2*phi))*cos(theta) - l2*(-cos(2*phi))*sin(theta))*ft;

        (l1*(-cos(2*phi))*sin(theta) + l2*(-sin(2*phi))*cos(theta))*fn + ...
        (l1*(sin(2*phi))*sin(theta) + l2*(-cos(2*phi))*cos(theta))*ft;

        l3*(-radius*sin(2*phi)*(-cos(2*phi)) + (len/2+radius*cos(2*phi))*(-sin(2*phi)))*fn + ...
        l3*(-radius*sin(2*phi)*(sin(2*phi)) + (len/2+radius*cos(2*phi))*(-cos(2*phi)))*ft;

        phi_dot ];

    dx_capsule2 = [
        (l1*0*cos(theta) - l2*(-1)*sin(theta))*fn + (l1*1*cos(theta) - l2*0*sin(theta))*ft;
        (l1*0*sin(theta) + l2*(-1)*cos(theta))*fn + (l1*1*sin(theta) + l2*0*cos(theta))*ft;
        l3*(-radius*0 + radius*(cos(phi)/sin(phi))*(-1))*fn + l3*(-radius*1 + radius*(cos(phi)/sin(phi))*0)*ft;
        phi_dot ];

    dx_capsule3 = [
        (l1*(cos(2*phi))*cos(theta) - l2*(sin(2*phi))*sin(theta))*fn + ...
        (l1*(-sin(2*phi))*cos(theta) - l2*(cos(2*phi))*sin(theta))*ft;
        (l1*(cos(2*phi))*sin(theta) + l2*(sin(2*phi))*cos(theta))*fn + ...
        (l1*(-sin(2*phi))*sin(theta) + l2*(cos(2*phi))*cos(theta))*ft;
        l3*(radius*sin(2*phi)*(cos(2*phi)) - (len/2+radius*cos(2*phi))*(sin(2*phi)))*fn + ...
        l3*(radius*sin(2*phi)*(-sin(2*phi)) - (len/2+radius*cos(2*phi))*(cos(2*phi)))*ft;
        phi_dot ];

    dx_capsule4 = [
        (l1*0*cos(theta) - l2*1*sin(theta))*fn + (l1*(-1)*cos(theta) - l2*0*sin(theta))*ft;
        (l1*0*sin(theta) + l2*1*cos(theta))*fn + (l1*(-1)*sin(theta) + l2*0*cos(theta))*ft;
        l3*(radius*0 - radius*(cos(phi)/sin(phi))*1)*fn + l3*(radius*(-1) - radius*(cos(phi)/sin(phi))*0)*ft;
        phi_dot ];

    expr_capsule = if_else(phi >= 0 & phi < pi/4 | phi >= 7*pi/4 & phi < 2*pi, dx_capsule1, ...
                     if_else(phi >= pi/4 & phi < 3*pi/4, dx_capsule2, ...
                     if_else(phi >= 3*pi/4 & phi < 5*pi/4, dx_capsule3, dx_capsule4)));

    % === 3) Prism piecewise cases ===
    dx_prism1 = [
        (l1*cos(theta)*(-1) - l2*sin(theta)*0)*fn + (l1*cos(theta)*0 - l2*sin(theta)*(-1))*ft;
        (l1*sin(theta)*(-1) + l2*cos(theta)*0)*fn + (l1*sin(theta)*0 + l2*cos(theta)*(-1))*ft;
        l3*(-len/2*tan(phi)*(-1) + len/2*0)*fn + l3*(-len/2*tan(phi)*0 + len/2*(-1))*ft;
        phi_dot ];

    dx_prism2 = [
        (l1*cos(theta)*0 - l2*sin(theta)*(-1))*fn + (l1*cos(theta)*1 - l2*sin(theta)*0)*ft;
        (l1*sin(theta)*0 + l2*cos(theta)*(-1))*fn + (l1*sin(theta)*1 + l2*cos(theta)*0)*ft;
        l3*((-wid/2)*0 + wid/2*(cos(phi)/sin(phi)*(-1)))*fn + l3*((-wid/2)*1 + wid/2*(cos(phi)/sin(phi)*0))*ft;
        phi_dot ];

    dx_prism3 = [
        (l1*cos(theta)*1 - l2*sin(theta)*0)*fn + (l1*cos(theta)*0 - l2*sin(theta)*1)*ft;
        (l1*sin(theta)*1 + l2*cos(theta)*0)*fn + (l1*sin(theta)*0 + l2*cos(theta)*1)*ft;
        l3*(len/2*tan(phi)*1 + (-len/2*0))*fn + l3*(len/2*tan(phi)*0 + (-len/2*1))*ft;
        phi_dot ];

    dx_prism4 = [
        (l1*cos(theta)*0 - l2*sin(theta)*1)*fn + (l1*(-1)*cos(theta) - l2*sin(theta)*0)*ft;
        (l1*sin(theta)*0 + l2*cos(theta)*1)*fn + (l1*(-1)*sin(theta) + l2*cos(theta)*0)*ft;
        l3*(wid/2*0 + (-wid/2*(cos(phi)/sin(phi)*1)))*fn + l3*(wid/2*(-1) + (-wid/2*(cos(phi)/sin(phi)*0)))*ft;
        phi_dot ];

    expr_prism = if_else(phi >= 0 & phi < alpha | phi >= 2*pi - alpha & phi < 2*pi, dx_prism1, ...
                  if_else(phi >= alpha & phi < pi - alpha, dx_prism2, ...
                  if_else(phi >= pi - alpha & phi < pi + alpha, dx_prism3, dx_prism4)));

    % === 4) Use MATLAB if/elseif ===
    if object_shape == "rectangular_capsule_prism"
        rhs = expr_capsule;
    elseif object_shape == "rectangular_prism"
        rhs = expr_prism;
    else
        error('Unknown object shape!');
    end

    % === 5) Jacobians ===
    A_sym = jacobian(rhs, x_vec);
    B_sym = jacobian(rhs, u_vec);

    A_fun = Function('A_fun', {x_vec, u_vec}, {A_sym});
    B_fun = Function('B_fun', {x_vec, u_vec}, {B_sym});

end
