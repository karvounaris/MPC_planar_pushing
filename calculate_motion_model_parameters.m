%========================================================================%
% This function calculates necessary parameters for the motion equation, given:
%     - u (control input)
%     - theta
%     - len
%     - radius
%     - phi
%========================================================================%

function wrench = calculate_motion_model_parameters(u, theta, len, radius, phi)

% Set the control input
u_f = u(1:2);

% Calculate rotation matrix
R = [cos(theta) -sin(theta) 0;
     sin(theta) cos(theta) 0;
     0 0 1];

% Calculate wrench from all the contact points to the object relative to
% F_a
[x_c, y_c, ~, n_c, t_c] = calculate_r_c(phi, len, radius);

J_c =  [1 0 -y_c;
        0 1 x_c];

N = J_c' * n_c;
T = J_c' * t_c;

B = [N T];

wrench = R * B * u_f;

end






















% %========================================================================%
% % This function calculates necessary parameters for the motion equation, given:
% %     - u (control input)
% %     - theta
% %     - len
% %     - radius
% %     - C (number of contact points)
% %     - contact_points (structure of contact points)
% %     - mass
% %     - N (auxiliary matrix)
% %     - T (auxiliary matrix)
% %     - mass
% %     - mu_ground (coefficient of friction)
% %     - contact_area (are of the object in touch with the ground)
% %     - dp (slight change of the position and orientation)
% %     - phi_C_start (first value of phi angle)
% %========================================================================%
% 
% function [R, contact_points, w, ground_friction] = calculate_motion_model_parameters(u, theta, len, radius, C, ...
%                                                     contact_points, N, T, mass, mu_ground, contact_area, dp, phi)
% 
% % Set the control input
% u_f = u(1:2*C);
% 
% % Calculate rotation matrix
% R = [cos(theta) -sin(theta) 0;
%      sin(theta) cos(theta) 0;
%      0 0 1];
% 
% % Calculate wrench from all the contact points to the object relative to
% % F_a
% for i = 1:C
%     [x_C, y_C, ~, n_C, t_C] = calculate_r_c(phi(i), len, radius);
%     contact_points{i,1} = x_C;
%     contact_points{i,2} = y_C;
%     contact_points{i,3} = phi(i);
%     contact_points{i,4} = [1 0 -y_C;
%                            0 1 x_C];
%     contact_points{i,5} = n_C;
%     contact_points{i,6} = t_C;
% end
% 
% for i = 1:length(N(1,:))
%     N(:,i) = contact_points{i,4}' * contact_points{i,5};
%     T(:,i) = contact_points{i,4}' * contact_points{i,6};
% end
% B = [N,T];
% 
% w = B * u_f;
% w = R * w;
% 
% % Calculate ground friction
% ground_friction = calculate_friction_with_ground(mass, mu_ground, contact_area, dp);
% 
% end