%=======================================================================%
% This function provides the shape of the capsule object for plotting, given:
%     - len    (length of the rectagnular side of the object)
%     - radius    (radius of semi-circular sides of the object)
%     - x   (x coordinate of center of mass)
%     - y   (y coordinate of center of mass)
%     - theta   (angle of orientation)
%=======================================================================%

function vertices = get_capsule_shape(len, radius, x, y, theta)
    half_len = len / 2;
   
    rect_x_1 = [half_len, -half_len];
    rect_x_2 = [-half_len, half_len];
    rect_y_1 = [radius, radius];
    rect_y_2 = [-radius, -radius];
    
    phi_circle = linspace(-pi/2, pi/2, 100);
    circle_x1 = half_len + radius * cos(phi_circle);
    circle_y1 = radius * sin(phi_circle);
    
    phi_circle = linspace(pi/2, 3* pi/2, 100);
    circle_x2 = -half_len + radius * cos(phi_circle);
    circle_y2 = radius * sin(phi_circle);
    
    capsule_x = [circle_x1, rect_x_1, circle_x2, rect_x_2];
    capsule_y = [circle_y1, rect_y_1, circle_y2, rect_y_2];
    
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    rotated_vertices = R * [capsule_x; capsule_y];
    
    vertices = rotated_vertices' + [x, y];
end