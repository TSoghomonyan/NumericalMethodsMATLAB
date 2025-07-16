% Symmetric system Ax = b solved via modified LU / LDLᵗ-like factorization
n = 3;
A = [9, -3, 1; -3, 2, 0; 1, 0, 4];
b = [5; 3; 2];

% Check symmetry and leading minors
isValid = true;
for i = 1:n
    if det(A(1:i,1:i)) == 0
        disp("Submatrix A(1:i,1:i) is singular");
        isValid = false;
    end
end

if ~issymmetric(A) || ~isValid
    disp("Matrix A is not symmetric or not positive definite.");
else
    % LU decomposition (Doolittle style)
    L = zeros(n);
    U = eye(n);

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

    disp("L:");
    disp(L);
    disp("U:");
    disp(U);

    % Extract diagonal D from L
    D = eye(n);
    for i = 1:n
        D(i,i) = L(i,i);
        L(i,i) = 1;  % Normalize L
    end

    % Normalize lower triangle of L by dividing by D(j,j)
    for i = 1:n
        for j = 1:i-1
            L(i,j) = L(i,j) / D(j,j);
        end
    end

    disp("Normalized L:");
    disp(L);
    disp("Diagonal D:");
    disp(D);

    % Take square root of D diagonals (Cholesky-like step)
    for i = 1:n
        D(i,i) = sqrt(D(i,i));
    end

    disp("sqrt(D):");
    disp(D);

    % Form Y = D * U
    Y = D * U;

    % Solve Yᵗ y = b using forward substitution
    S = Y';
    y = zeros(n,1);
    for i = 1:n
        s = 0;
        for j = 1:i-1
            s = s + S(i,j) * y(j);
        end
        y(i) = (b(i) - s) / S(i,i);
    end

    disp("Intermediate y:");
    disp(y);

    % Solve Y x = y using backward substitution
    x = zeros(n,1);
    for i = n:-1:1
        s = 0;
        for j = i+1:n
            s = s + Y(i,j) * x(j);
        end
        x(i) = (y(i) - s) / Y(i,i);
    end

    disp("Solution x:");
    disp(x);
end
