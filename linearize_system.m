function [A, B] = linearize_system(x_star, u_star, L, radius, len, wid)
    % Define symbolic variables for state and input
    alpha = atan(wid/len);
    syms x y theta phi fn ft phi_dot L l1 l2 l3 len wid
    x_vec = [x; y; theta; phi];
    u_vec = [fn; ft; phi_dot];

    L = [l1 0 0;
         0 l2 0;
         0 0 l3];
    
    % Define symbolic expressions for each case in the dynamics capsule
    % dx_case1 = [
    %     (L(1,1)*(-cos(2*phi))*cos(theta) - L(2,2)*(-sin(2*phi))*sin(theta))*fn + ...
    %     (L(1,1)*(sin(2*phi))*cos(theta) - L(2,2)*(-cos(2*phi))*sin(theta))*ft;
    % 
    %     (L(1,1)*(-cos(2*phi))*sin(theta) + L(2,2)*(-sin(2*phi))*cos(theta))*fn + ...
    %     (L(1,1)*(sin(2*phi))*sin(theta) + L(2,2)*(-cos(2*phi))*cos(theta))*ft;
    % 
    %     L(3,3)*(-radius*sin(2*phi)*(-cos(2*phi)) + (len/2+radius*cos(2*phi))*(-sin(2*phi)))*fn + ...
    %     L(3,3)*(-radius*sin(2*phi)*(sin(2*phi)) + (len/2+radius*cos(2*phi))*(-cos(2*phi)))*ft;
    % 
    %     phi_dot];
    % 
    % dx_case2 = [
    %     (L(1,1)*(0)*cos(theta) - L(2,2)*(-1)*sin(theta))*fn + ...
    %     (L(1,1)*(1)*cos(theta) - L(2,2)*(0)*sin(theta))*ft;
    % 
    %     (L(1,1)*(0)*sin(theta) + L(2,2)*(-1)*cos(theta))*fn + ...
    %     (L(1,1)*(1)*sin(theta) + L(2,2)*(0)*cos(theta))*ft;
    % 
    %     L(3,3)*(-radius*(0) + radius*(cos(phi)/sin(phi))*(-1))*fn + ...
    %     L(3,3)*(-radius*(1) + radius*(cos(phi)/sin(phi))*(0))*ft;
    % 
    %     phi_dot];
    % 
    % dx_case3 = [
    %     (L(1,1)*(cos(2*phi))*cos(theta) - L(2,2)*(sin(2*phi))*sin(theta))*fn + ...
    %     (L(1,1)*(-sin(2*phi))*cos(theta) - L(2,2)*(cos(2*phi))*sin(theta))*ft;
    % 
    %     (L(1,1)*(cos(2*phi))*sin(theta) + L(2,2)*(sin(2*phi))*cos(theta))*fn + ...
    %     (L(1,1)*(-sin(2*phi))*sin(theta) + L(2,2)*(cos(2*phi))*cos(theta))*ft;
    % 
    %     L(3,3)*(radius*sin(2*phi)*(cos(2*phi)) - (len/2+radius*cos(2*phi))*(sin(2*phi)))*fn + ...
    %     L(3,3)*(radius*sin(2*phi)*(-sin(2*phi)) - (len/2+radius*cos(2*phi))*(cos(2*phi)))*ft;
    % 
    %     phi_dot];
    % 
    % dx_case4 = [
    %     (L(1,1)*(0)*cos(theta) - L(2,2)*(1)*sin(theta))*fn + ...
    %     (L(1,1)*(-1)*cos(theta) - L(2,2)*(0)*sin(theta))*ft;
    % 
    %     (L(1,1)*(0)*sin(theta) + L(2,2)*(1)*cos(theta))*fn + ...
    %     (L(1,1)*(-1)*sin(theta) + L(2,2)*(0)*cos(theta))*ft;
    %
    %     L(3,3)*(radius*(0) - radius*(cos(phi)/sin(phi))*(1))*fn + ...
    %     L(3,3)*(radius*(-1) - radius*(cos(phi)/sin(phi))*(0))*ft;
    % 
    %     phi_dot];

    % % Select the appropriate case based on phi value
    % phi_val = x_star(4);
    % if (phi_val >= 0 && phi_val < pi/4) || (phi_val >= 7*pi/4 && phi_val < 2*pi)
    %     f_sym = dx_case1;
    % elseif phi_val >= pi/4 && phi_val < 3*pi/4
    %     f_sym = dx_case2;
    % elseif phi_val >= 3*pi/4 && phi_val < 5*pi/4
    %     f_sym = dx_case3;
    % elseif phi_val >= 5*pi/4 && phi_val < 7*pi/4
    %     f_sym = dx_case4;
    % end

    x_c = len/2;
    y_c = (len/2)*tan(phi);
    n_c = [-1; 0];
    t_c = [0; -1];
    % Define symbolic expressions for each case in the dynamics square
    dx_case1 = [(L(1,1)*cos(theta)*(-1) - L(2,2)*sin(theta)*0)*fn + ...
                (L(1,1)*cos(theta)*0 - L(2,2)*sin(theta)*(-1))*ft; ...
        
                (L(1,1)*sin(theta)*(-1) + L(2,2)*cos(theta)*0)*fn + ...
                (L(1,1)*sin(theta)*0 + L(2,2)*cos(theta)*(-1))*ft; ...
                
                (L(3,3)*(-len/2*tan(phi)*(-1) + (len/2*0)))*fn + ...
                (L(3,3)*(-len/2*tan(phi)*0 + (len/2*(-1))))*ft;
                
                phi_dot];

    dx_case2 = [(L(1,1)*cos(theta)*0 - L(2,2)*sin(theta)*(-1))*fn + ...
                (L(1,1)*cos(theta)*1 - L(2,2)*sin(theta)*0)*ft; ...
        
                (L(1,1)*sin(theta)*0 + L(2,2)*cos(theta)*(-1))*fn + ...
                (L(1,1)*sin(theta)*1 + L(2,2)*cos(theta)*0)*ft; ...
                
                (L(3,3)*((-wid/2)*0 + (wid/2*(cos(phi)/sin(phi)*(-1)))))*fn + ...
                (L(3,3)*((-wid/2)*1 + (wid/2*(cos(phi)/sin(phi)*(0)))))*ft;
                
                phi_dot];

    dx_case3 = [(L(1,1)*cos(theta)*1 - L(2,2)*sin(theta)*0)*fn + ...
                (L(1,1)*cos(theta)*0 - L(2,2)*sin(theta)*1)*ft; ...
        
                (L(1,1)*sin(theta)*1 + L(2,2)*cos(theta)*0)*fn + ...
                (L(1,1)*sin(theta)*0 + L(2,2)*cos(theta)*1)*ft; ...
                
                (L(3,3)*(len/2*tan(phi)*1 + (-len/2*0)))*fn + ...
                (L(3,3)*(len/2*tan(phi)*0 + (-len/2*1)))*ft;
                
                phi_dot];

    dx_case4 = [(L(1,1)*cos(theta)*0 - L(2,2)*sin(theta)*1)*fn + ...
                (L(1,1)*cos(theta)*(-1) - L(2,2)*sin(theta)*0)*ft; ...
        
                (L(1,1)*sin(theta)*0 + L(2,2)*cos(theta)*1)*fn + ...
                (L(1,1)*sin(theta)*(-1) + L(2,2)*cos(theta)*0)*ft; ...
                
                (L(3,3)*((wid/2)*0 + (-wid/2*(cos(phi)/sin(phi)*1))))*fn + ...
                (L(3,3)*((wid/2)*(-1) + (-wid/2*(cos(phi)/sin(phi)*0))))*ft;
                
                phi_dot];

    phi_val = x_star(4);
    
    if ((phi_val >= 0) && (phi_val < alpha)) || ((phi_val >= 2*pi - alpha) && (phi_val < 2*pi))
        f_sym = dx_case1;
    elseif (phi_val >= alpha) && (phi_val < pi - alpha)
        f_sym = dx_case2;
    elseif (phi_val >= pi - alpha) && (phi_val < pi + alpha)
        f_sym = dx_case3;
    elseif (phi_val >= pi + alpha) && (phi_val < 2*pi - alpha)
        f_sym = dx_case4;
    end
    
    % Calculate the Jacobians of f_sym with respect to x and u
    A_sym = jacobian(f_sym, x_vec)
    B_sym = jacobian(f_sym, u_vec)

    % Substitute x_star and u_star into the Jacobian matrices
    % A = double(subs(A_sym, {x, y, theta, phi, fn, ft, phi_dot}, ...
    %                 {x_star(1), x_star(2), x_star(3), x_star(4), u_star(1), u_star(2), u_star(3)}));
    % B = double(subs(B_sym, {x, y, theta, phi, fn, ft, phi_dot}, ...
    %                 {x_star(1), x_star(2), x_star(3), x_star(4), u_star(1), u_star(2), u_star(3)}));
end
