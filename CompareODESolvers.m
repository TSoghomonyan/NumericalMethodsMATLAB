%===========================================================
% Title  : Solving a Cauchy Problem Using Euler and RK4 Methods
% 
% Problem:
%   Solve the ODE y' = -x*y, y(0) = 1 over the interval [0, 3]
%   using:
%     - Exact solution
%     - Explicit Euler method
%     - Implicit Euler method
%     - Runge-Kutta method of 4th order (RK4)
% 
% Exact solution: y(x) = exp(-x^2 / 2)
%===========================================================

% Exact solution
x_exact = linspace(0, 3, 30);
y_exact = exp(-x_exact.^2 / 2);
plot(x_exact, y_exact, 'b', 'LineWidth', 1.5);
hold on;

% Step size and interval
n = 30;
h = 3 / n;
x = 0:h:3;

%-----------------------------------------------------------
% Explicit Euler Method
%-----------------------------------------------------------
Y_explicit = zeros(1, length(x));
Y_explicit(1) = 1;

for i = 2:length(x)
    Y_explicit(i) = Y_explicit(i-1) + h * (-x(i-1) * Y_explicit(i-1));
end

plot(x, Y_explicit, 'r--', 'LineWidth', 1.5);

%-----------------------------------------------------------
% Implicit Euler Method
%-----------------------------------------------------------
Y_implicit = zeros(1, length(x));
Y_implicit(1) = 1;

for i = 2:length(x)
    Y_implicit(i) = Y_implicit(i-1) / (1 + h * x(i));
end

plot(x, Y_implicit, 'g--', 'LineWidth', 1.5);

%-----------------------------------------------------------
% Runge-Kutta 4th Order Method (RK4)
%-----------------------------------------------------------
Y_rk4 = zeros(1, length(x));
Y_rk4(1) = 1;

for i = 2:length(x)
    k1 = -x(i-1) * Y_rk4(i-1);
    k2 = -(x(i-1) + h/2) * (Y_rk4(i-1) + h/2 * k1);
    k3 = -(x(i-1) + h/2) * (Y_rk4(i-1) + h/2 * k2);
    k4 = -(x(i-1) + h) * (Y_rk4(i-1) + h * k3);
    
    Y_rk4(i) = Y_rk4(i-1) + (h / 6) * (k1 + 2*k2 + 2*k3 + k4);
end

plot(x, Y_rk4, 'k-.', 'LineWidth', 1.5);

%-----------------------------------------------------------
% Graph formatting
%-----------------------------------------------------------
legend('Exact Solution', ...
       'Explicit Euler', ...
       'Implicit Euler', ...
       'Runge-Kutta 4th Order');
xlabel('x');
ylabel('y(x)');
title('Solving y'' = -x*y using Euler and RK4 Methods');
grid on;
hold off;
