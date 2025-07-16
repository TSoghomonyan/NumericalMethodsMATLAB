% LU decomposition without pivoting using elementary matrices
% First: fixed 3x3 matrix for educational example

A = [3, 18, -6; -5, -23, 31; 9, 65, 13];
U = A;
n = size(A,1);
L = eye(n);
CurrentL = L;

for i = 1:n
    CurrentL = eye(n);
    for j = i:n
        if j == i
            CurrentL(i,j) = CurrentL(i,j) / U(i,i);
        else
            CurrentL(j,i) = -U(j,i) / U(i,i);
        end
    end
    U = CurrentL * U;
    L = CurrentL * L;
end

L = inv(L);  % Reversing the product of elementary matrices
disp("L matrix:");
disp(L);
disp("U matrix:");
disp(U);

%% Now generalize to random nxn matrix

n = 5;
isInvertible = false;

while ~isInvertible
    A = round(rand(n,n) * 10);
    isInvertible = true;

    % Check leading principal minors
    for i = 1:n
        if det(A(1:i,1:i)) == 0
            disp("Submatrix is singular, retrying...");
            isInvertible = false;
            break;
        end
    end
end

U = A;
L = eye(n);

for i = 1:n
    CurrentL = eye(n);
    for j = i:n
        if j == i
            CurrentL(i,j) = CurrentL(i,j) / U(i,i);
        else
            CurrentL(j,i) = -U(j,i) / U(i,i);
        end
    end
    U = CurrentL * U;
    L = CurrentL * L;
end

L = inv(L);

disp("LU decomposition successful.");
disp("Matrix A:");
disp(A);
disp("Lower triangular matrix L:");
disp(L);
disp("Upper triangular matrix U:");
disp(U);

% Optional: Check if A ≈ L*U
disp("L*U:");
disp(L*U);
