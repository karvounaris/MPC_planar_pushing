function [objGrad, cGrad, ceqGrad, hessian] = compute_derivatives(z, ~, lambda, data)
% Computes analytical derivatives for:
% - objective gradient
% - nonlinear constraint Jacobian
% - Hessian of Lagrangian (placeholder here)
%
% Inputs:
%   z      - decision variable vector
%   lambda - Lagrange multipliers for constraints (struct with fields .ineq and .eq)
%   data   - struct with parameters (Q, QN, R, mu, L, dt, nx, nu, N, x_star, len, radius, wid, object_shape)
%
% Outputs:
%   objGrad - gradient of objective (nx*(N+1)+nu*N+N,1)
%   cGrad   - Jacobian of nonlinear inequalities (empty here)
%   ceqGrad - Jacobian of nonlinear equalities (nx*N + N, nz sparse)
%   hessian - Hessian of Lagrangian (nz x nz sparse) (empty here)

Q = data.Q; QN = data.QN; R = data.R; w_eps0 = data.w_eps0; k_eps = data.k_eps;
mu = data.mu; Lmat = data.L; dt = data.dt; nx = data.nx; nu = data.nu; N = data.N;
x_star = data.x_star; len = data.len; radius = data.radius; wid = data.wid; object_shape = data.object_shape;

nX = nx*(N+1);
nU = nu*N;
nE = N;
nz = length(z);

X = reshape(z(1:nX), nx, []);
U = reshape(z(nX+1:nX+nU), nu, []);
Eps = z(nX+nU+1:end);

objGrad = zeros(nz,1);
cGrad = sparse(0, nz); % No nonlinear inequalities
ceqGrad = sparse(nx*N + N, nz);

% --- Gradient of objective ---

% States
for i=1:N
    dx = X(:,i) - x_star(:,i);
    objGrad((i-1)*nx + (1:nx)) = 2 * Q * dx;
end
dxN = X(:,N+1) - x_star(:,N+1);
objGrad(N*nx + (1:nx)) = 2 * QN * dxN;

% Controls
for i=1:N
    u_idx = nX + (i-1)*nu;
    objGrad(u_idx + (1:nu)) = 2 * R * U(:,i);
end

% Slack variables
for i=1:N
    eps_idx = nX + nU + i;
    w_eps_i = w_eps0 * exp(-k_eps^(i-1));
    objGrad(eps_idx) = 2 * w_eps_i * Eps(i);
end

% --- Jacobian of nonlinear equality constraints ---

% Constraint vector ordering:
% 1 to nx*N: dynamics constraints per step (each length nx)
% nx*N+1 to nx*N+N: complementarity constraints per step (scalar)

for i=1:N
    % Indices for variables
    x_i_idx = (i-1)*nx + (1:nx);
    x_ip1_idx = i*nx + (1:nx);
    u_i_idx = nX + (i-1)*nu + (1:nu);
    eps_i_idx = nX + nU + i;

    % Extract relevant states and inputs
    theta = X(3,i);
    phi = X(4,i);
    [x_c, y_c, ~, n_c, t_c] = calculate_r_c(phi, len, radius, wid, object_shape);

    Rot = [cos(theta), -sin(theta), 0;
           sin(theta),  cos(theta), 0;
           0,           0,          1];
    Jc = [1 0 -y_c; 0 1 x_c];
    B = [Jc' * n_c, Jc' * t_c];

    % Derivatives needed:

    % dRot/dtheta
    dRot_dtheta = [-sin(theta), -cos(theta), 0;
                    cos(theta), -sin(theta), 0;
                    0,          0,          0];

    % dJc/dphi (approximate as zeros for now)
    dyc_dphi = 0; % You can add analytical derivatives here if known
    dxc_dphi = 0;

    dJc_dphi = [0 0 -dyc_dphi; 0 0 dxc_dphi]; % zeros here

    % dB/dphi
    dB_dphi = [dJc_dphi' * n_c, Jc' * zeros(2,1)]; % approx zero vector

    % Compute partials of fx = [Rot * L * B * [fn; ft]; phidot+ - phidot-]

    % State derivatives:

    % dfx/dtheta (4x1)
    dfx_dtheta = zeros(nx,1);
    dfx_dtheta(1:3) = dRot_dtheta * Lmat * B * U(1:2,i);
    dfx_dtheta(4) = 0;

    % dfx/dphi (4x1)
    dfx_dphi = zeros(nx,1);
    dfx_dphi(1:3) = Rot * Lmat * dB_dphi * U(1:2,i);
    dfx_dphi(4) = 0;

    % Jacobian of fx wrt x_i (nx x nx)
    dfx_dx_i = zeros(nx,nx);
    dfx_dx_i(:,3) = dfx_dtheta;
    dfx_dx_i(:,4) = dfx_dphi;

    % Jacobian of fx wrt u_i (nx x nu)
    dfx_du_i = zeros(nx,nu);
    dfx_du_i(1:3,1:2) = Rot * Lmat * B;
    dfx_du_i(4,3) = 1;
    dfx_du_i(4,4) = -1;

    % Dynamics constraints: ceq_dyn = x_{i+1} - x_i - dt*fx = 0
    ceqGrad((i-1)*nx + (1:nx), x_i_idx) = -eye(nx) - dt * dfx_dx_i;
    ceqGrad((i-1)*nx + (1:nx), x_ip1_idx) = eye(nx);
    ceqGrad((i-1)*nx + (1:nx), u_i_idx) = -dt * dfx_du_i;

    % Complementarity constraint:
    fn = U(1,i);
    ft = U(2,i);
    phi_dot_plus = U(3,i);
    phi_dot_minus = U(4,i);

    lambda_minus = mu*fn + ft;
    lambda_plus = mu*fn - ft;

    % Partial derivatives of complementarity constraint:
    % ceq_comp = lambda_minus * phi_dot_plus + lambda_plus * phi_dot_minus + Eps(i) = 0

    dcomp_dfn = mu * phi_dot_plus + mu * phi_dot_minus;
    dcomp_dft = phi_dot_plus - phi_dot_minus;
    dcomp_dphi_dot_plus = lambda_minus;
    dcomp_dphi_dot_minus = lambda_plus;
    dcomp_deps = 1;

    ceqGrad(nx*N + i, u_i_idx(1)) = dcomp_dfn;
    ceqGrad(nx*N + i, u_i_idx(2)) = dcomp_dft;
    ceqGrad(nx*N + i, u_i_idx(3)) = dcomp_dphi_dot_plus;
    ceqGrad(nx*N + i, u_i_idx(4)) = dcomp_dphi_dot_minus;
    ceqGrad(nx*N + i, eps_i_idx) = dcomp_deps;
end

% Hessian (optional, empty to let Knitro approximate)
hessian = sparse(nz, nz);

end
