% Direct implementation of Gauss-Seidel method for solving Ax = b

A = [2, 1, 0; 1, 4, 2; -1, 3, 5];
b = [1; 1; 0];
n = length(b);

% Initial guess
x = zeros(n, 1);

% Number of iterations
num_iter = 6;

for z = 1:num_iter
    for i = 1:n
        s = 0;
        for j = 1:n
            if i == j
                continue
            end
            s = s + A(i,j) * x(j);
        end
        % Gauss-Seidel formula
        x(i) = (b(i) - s) / A(i,i);
    end
    disp(['Iteration ', num2str(z), ':']);
    disp(x);
end
