%=======================================================================%
% This function provides the shape of a circle for plotting, given:
%     - radius (radius of the circle)
%     - x   (x coordinate of center)
%     - y   (y coordinate of center)
%=======================================================================%

function vertices = get_circle_shape(radius, x, y)
    theta = linspace(0, 2*pi, 100);
    circle_x = radius * cos(theta);
    circle_y = radius * sin(theta);
    vertices = [circle_x' + x, circle_y' + y];
end
