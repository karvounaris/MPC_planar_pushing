%=======================================================================%
% This function calculates the inertia matrix of a rectangular capsule 
% prism given the:
%     - length
%     - width
%     - height
%     - radius
%     - rectangular_prism_mass
%     - half_cylinder_mass
%     - object_shape ("rectangular_capsule_prism" or "rectangular_prism")
%=======================================================================%

function inertia_matrix = calculate_inertia_matrix(length, width, height, ...
                              radius, rectangular_prism_mass, half_cylinder_mass, mass, ...
                              object_shape)
    
    if object_shape == "rectangular_capsule_prism"
        I_xx = (1/12) * rectangular_prism_mass * (height^2 + width^2) + ...
               (1/2) * 2 * half_cylinder_mass * (3*radius^2 + height^2);
        I_yy = (1/12) * rectangular_prism_mass * (height^2 + length^2) + ...
               (1/2) * 2 * half_cylinder_mass * (3*radius^2 + height^2) + ...
               2 * half_cylinder_mass * (length/2)^2;
        I_zz = (1/12) * rectangular_prism_mass * (width^2 + length^2) + ...
               (1/2) * 2 * half_cylinder_mass * radius^2 + ...
               2 * half_cylinder_mass * (length/2)^2;
    elseif object_shape == "rectangular_prism"
        I_xx = (1/12) * mass * (height^2 + width^2);
        I_yy = (1/12) * mass * (height^2 + length^2);
        I_zz = (1/12) * mass * (width^2 + length^2);
    else
        disp("Invalid object shape");
    end
    
    inertia_matrix = [I_xx 0 0;
                      0 I_yy 0;
                      0 0 I_zz];
end