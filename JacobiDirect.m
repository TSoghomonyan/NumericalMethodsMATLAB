% Direct implementation of the Jacobi method for solving Ax = b

A = [2, 1, 0; 1, 4, 2; -1, 3, 5];
b = [1; 1; 0];
n = 3;

% Initial guess
x = zeros(n,1);

% Number of iterations
num_iter = 6;

for z = 1:num_iter
    xNew = zeros(n,1);
    for i = 1:n
        s = 0;
        for j = 1:n
            if i == j
                continue
            end
            s = s + A(i,j) * x(j);
        end
        xNew(i) = (b(i) - s) / A(i,i);
    end
    disp(['Iteration ', num2str(z), ':']);
    disp(x);
    x = xNew;
end
