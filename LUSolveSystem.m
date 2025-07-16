% Solve Ax = b using LU decomposition (Doolittles method)
% Then forward and backward substitution

disp("Example 1: 3x3 system");

n = 3;
A = [3, 18, -6; -5, -23, 31; 9, 65, 13];
b = [9; 7; 2];

L = zeros(n);
U = eye(n);

% LU decomposition
for i = 1:n
    for j = 1:n
        s = A(i,j);
        if i >= j
            for k = 1:j-1
                s = s - L(i,k) * U(k,j);
            end
            L(i,j) = s / U(j,j);
        else
            for k = 1:i-1
                s = s - L(i,k) * U(k,j);
            end
            U(i,j) = s / L(i,i);
        end
    end
end

disp("L matrix:");
disp(L);
disp("U matrix:");
disp(U);

% Forward substitution: Lz = b
z = zeros(n,1);
for i = 1:n
    s = 0;
    for j = 1:i-1
        s = s + L(i,j) * z(j);
    end
    z(i) = (b(i) - s) / L(i,i);
end

disp("Intermediate vector z:");
disp(z);

% Backward substitution: Ux = z
x = zeros(n,1);
for i = n:-1:1
    s = 0;
    for j = i+1:n
        s = s + U(i,j) * x(j);
    end
    x(i) = (z(i) - s) / U(i,i);
end

disp("Solution vector x:");
disp(x);

%% ========== Second Example ==========

disp("Example 2: 4x4 system");

n = 4;
A = [1, 5, 3, 4; 1, 8, 15, 16; 1, 12, 33, 42; 1, 9, 24, 48];
b = [29; 95; 217; 211];

L = zeros(n);
U = eye(n);

% LU decomposition
for i = 1:n
    for j = 1:n
        s = A(i,j);
        if i >= j
            for k = 1:j-1
                s = s - L(i,k) * U(k,j);
            end
            L(i,j) = s / U(j,j);
        else
            for k = 1:i-1
                s = s - L(i,k) * U(k,j);
            end
            U(i,j) = s / L(i,i);
        end
    end
end

disp("L matrix:");
disp(L);
disp("U matrix:");
disp(U);

% Forward substitution: Lz = b
z = zeros(n,1);
for i = 1:n
    s = 0;
    for j = 1:i-1
        s = s + L(i,j) * z(j);
    end
    z(i) = (b(i) - s) / L(i,i);
end

disp("Intermediate vector z:");
disp(z);

% Backward substitution: Ux = z
x = zeros(n,1);
for i = n:-1:1
    s = 0;
    for j = i+1:n
        s = s + U(i,j) * x(j);
    end
    x(i) = (z(i) - s) / U(i,i);
end

disp("Solution vector x:");
disp(x);
