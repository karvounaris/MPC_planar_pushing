%=======================================================================%
% This function provides the shape of a rectangle for plotting, given:
%     - len      (rectangle length along its local x-axis)
%     - wid      (rectangle width along its local y-axis)
%     - x        (x coordinate of center of mass)
%     - y        (y coordinate of center of mass)
%     - theta    (angle of orientation, in radians)
%=======================================================================%
function vertices = get_rectangle_shape(len, wid, x, y, theta)
    % Half-dimensions
    half_len = len / 2;
    half_wid = wid / 2;

    % Rectangle corners in local coordinates (centered at origin)
    corners_x = [ half_len,  half_len, -half_len, -half_len];
    corners_y = [ half_wid, -half_wid, -half_wid,  half_wid];

    % Rotation matrix
    R = [cos(theta), -sin(theta);
         sin(theta),  cos(theta)];

    % Rotate and translate
    rotated_vertices = R * [corners_x; corners_y];
    vertices = rotated_vertices' + [x, y];
end
