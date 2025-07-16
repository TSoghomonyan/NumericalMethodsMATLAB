% Iterative solution of Ax = b using the Jacobi method

n = 3;
A = [3, 6, -9; 2, 2, 8; -9/8, -3, 1];
b = [1; 2; 9];

% Extract diagonal matrix D and compute its inverse
D = zeros(n, n);
for i = 1:n
    D(i,i) = 1 / A(i,i);
    A(i,i) = 0; % Remove diagonal from A, leaving L + U
end

% Compute iteration matrix B = -D^{-1}*(L + U)
B = -D * A;

% Display the iteration matrix
disp("Jacobi iteration matrix B = -D^{-1}(L+U):");
disp(B);

% Check if matrix is suitable for convergence
if IsCorrect(B)
    x = zeros(n,1);      % Initial guess
    for i = 1:400
        x = B * x + D * b;
    end
    disp("Approximate solution after 400 iterations:");
    disp(x);
else
    disp("Matrix A is not suitable for iterative methods.");
end
