%==========================================================
% Title: Numerical solution of 1D wave equation using finite differences
%
% Problem:
%   Solve u_tt = a^2 * u_xx + f(x)
%   with initial and boundary conditions:
%       u(x,0) = sin(pi * x)
%       u_t(x,0) = 0
%       u(0,t) = 0, u(L,t) = 0
% Numerical method:
%   Explicit finite difference scheme with central differences in space and time
%==========================================================

% Parameters
M = 100;       % Number of spatial steps
N = 200;       % Number of time steps
L = 1;         % Length of spatial domain
a = 1;         % Wave speed
T = 2 * L / a; % Total time
h = L / M;     % Spatial step size
tau = T / N;   % Time step size
r = a * tau / h;

% Functions
f = @(x) 0 .* x;                 % Source term f(x)
phi = @(x) sin(pi .* x);         % Initial displacement u(x,0)
phiSecondGrade = @(x) -(pi^2) * sin(pi .* x); % Second derivative of phi(x)
ksi = @(x) 0 .* x;               % Initial velocity u_t(x,0)
alpha = @(t) 0 .* t;             % Boundary condition at x=0
beta = @(t) 0 .* t;              % Boundary condition at x=L

% Create grids
V = zeros(N+1, M+1);             % Solution matrix V(time, space)
t = linspace(0, T, N+1);
x = linspace(0, L, M+1);

% Initial condition at t=0
for j = 1:M+1
    V(1,j) = phi(x(j));
end

% Approximate solution at t = tau (second time layer)
for j = 2:M
    V(2,j) = tau * ksi(x(j)) + (tau^2 / 2) * (a^2 * phiSecondGrade(x(j)) + f(x(j))) + V(1,j);
end

% Boundary conditions for all time layers
for i = 2:N+1
    V(i,1) = alpha(t(i));
    V(i,M+1) = beta(t(i));
end

% Main loop - update interior points
for i = 2:N
    for j = 2:M
        V(i+1,j) = r^2 * (V(i,j+1) - 2*V(i,j) + V(i,j-1)) + 2*V(i,j) - V(i-1,j) + tau^2 * f(x(j));
    end
end

% Visualization of numerical solution
[X, T_grid] = meshgrid(x, t);
figure;
surf(X, T_grid, V);
title('Numerical Solution of Wave Equation');
xlabel('x');
ylabel('t');
zlabel('u(x,t)');
shading interp;

% Exact solution computation
U_exact = zeros(N+1, M+1);
for i = 1:N+1
    for j = 1:M+1
        U_exact(i,j) = sin(pi * x(j)) * cos(pi * t(i));
    end
end

% Visualization of exact solution
figure;
surf(X, T_grid, U_exact);
title('Exact Solution of Wave Equation');
xlabel('x');
ylabel('t');
zlabel('u(x,t)');
shading interp;
