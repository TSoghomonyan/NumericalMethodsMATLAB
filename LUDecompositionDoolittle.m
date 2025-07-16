% LU decomposition using Doolittles method (no pivoting)

disp("This example is from the lesson");
n = 3;
A = [3,18,-6; -5,-23,31; 9,65,13];

L = zeros(n);
U = zeros(n);

for i = 1:n
    % Compute U(i,j) for j = i to n
    for j = i:n
        s = 0;
        for k = 1:i-1
            s = s + L(i,k) * U(k,j);
        end
        U(i,j) = A(i,j) - s;
    end

    % Compute L(j,i) for j = i to n
    for j = i:n
        if i == j
            L(i,i) = 1;
        else
            s = 0;
            for k = 1:i-1
                s = s + L(j,k) * U(k,i);
            end
            L(j,i) = (A(j,i) - s) / U(i,i);
        end
    end
end

disp("Matrix A:");
disp(A);
disp("Lower triangular matrix L:");
disp(L);
disp("Upper triangular matrix U:");
disp(U);

%% General case for random n x n matrix

disp("Now a random n x n matrix");

n = 5;
A = round(rand(n,n) * 10);

L = zeros(n);
U = zeros(n);

for i = 1:n
    for j = i:n
        s = 0;
        for k = 1:i-1
            s = s + L(i,k) * U(k,j);
        end
        U(i,j) = A(i,j) - s;
    end

    for j = i:n
        if i == j
            L(i,i) = 1;
        else
            s = 0;
            for k = 1:i-1
                s = s + L(j,k) * U(k,i);
            end
            L(j,i) = (A(j,i) - s) / U(i,i);
        end
    end
end

disp("Random matrix A:");
disp(A);
disp("Computed L:");
disp(L);
disp("Computed U:");
disp(U);

% Optional: check reconstruction
disp("L*U:");
disp(L*U);
