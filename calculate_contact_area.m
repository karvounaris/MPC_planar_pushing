%=======================================================================%
% This function calculates the area of the object the is in touch with the
% ground given:
%     - length
%     - width
%     - radius
%     - object_shape ("rectangular_capsule_prism" or "rectangular_prism"
%                      or "circular_prism")
%=======================================================================%

function contact_area = calculate_contact_area(len, width, radius, ...
                                               object_shape)

if object_shape == "rectangular_capsule_prism"
    rectangular_contact_area = len * width;
    cylinder_contact_area = pi * radius^2;
elseif object_shape == "rectangular_prism"
    rectangular_contact_area = len * width;
    cylinder_contact_area = 0;
elseif object_shape == "circular_prism"
    rectangular_contact_area = 0;
    cylinder_contact_area = pi * radius^2;
else
    disp("Invalid object shape")
end

contact_area = rectangular_contact_area + cylinder_contact_area;
end