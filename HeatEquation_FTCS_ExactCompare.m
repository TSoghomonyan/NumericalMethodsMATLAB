%====================================================================
% Title  : Solving Heat Equation with Source Term Using FTCS Scheme
%
% Problem:
%   Solve the PDE: u_t = u_xx + x on the domain x ∈ [0, 1], t ∈ [0, 1]
%   with initial condition: u(x, 0) = sin(3πx/2)
%   and boundary condition derived from exact solution.
%
% Numerical Method:
%   Forward Time Centered Space (FTCS) scheme
%
% Exact Solution:
%   u(x,t) = x*t + exp(-(3π/2)^2 * t) * sin(3πx/2)
%====================================================================

% Grid parameters
M = 10;        % Number of spatial points
N = 200;       % Number of time points
L = 1;         % Length of spatial domain
T = 1;         % End time

h = L / (M - 1);
tau = T / (N - 1);
r = tau / h^2;

x = 0:h:L;
t = 0:tau:T;

% Stability check
if r > 0.5
    disp('Stability condition violated: r > 0.5');
end

% Initialize solution matrix
u = zeros(N, M);
u(1, :) = sin(3 * pi * x / 2);   % Initial condition

% Time-stepping using FTCS
for n = 1:N-1
    for m = 2:M-1
        u(n+1, m) = (1 - 2*r) * u(n, m) + ...
                    r * (u(n, m-1) + u(n, m+1)) + ...
                    tau * x(m);  % Source term
    end
    % Right boundary condition (derived from exact)
    u(n+1, M) = t(n+1) * h + u(n+1, M-1);
end

% Visualization of numerical solution
[X, T_grid] = meshgrid(x, t);
figure;
surf(X, T_grid, u);
title('Numerical Solution (FTCS)');
xlabel('x');
ylabel('t');
zlabel('u(x,t)');
shading interp;

% Exact solution
u_exact = X .* T_grid + exp(-(3*pi/2)^2 .* T_grid) .* sin(3*pi*X/2);

% Visualization of exact solution
figure;
surf(X, T_grid, u_exact);
title('Exact Solution');
xlabel('x');
ylabel('t');
zlabel('u(x,t)');
shading interp;
