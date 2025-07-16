% Gauss-Seidel Iteration for solving Ax = b

A = [2, 1, 0; 1, 4, 2; -1, 3, 5];
b = [1; 1; 0];
n = 3;

% Extract the strictly upper triangular part of A
U = zeros(n,n);
for i = 1:n
    for j = (i+1):n
        U(i,j) = A(i,j);
    end
end

% Compute L + D
LD = A - U;

% Compute inverse of (L + D)
LD_inv = inv(LD);

% Iteration matrix and vector
B = -LD_inv * U;
c = LD_inv * b;

% Display iteration matrix
disp("Iteration matrix B = -(L+D)^(-1)*U:");
disp(B);

% Check convergence
if IsCorrect(B)
    x = zeros(n,1);
    for i = 1:5
        x = B * x + c;
    end
    disp("Approximate solution after 5 iterations:");
    disp(x);
else
    disp("A's eigenvalues have absolute value >= 1; may not converge.");
end
