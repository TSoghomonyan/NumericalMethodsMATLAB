function ret = IsCorrect(A)
% IsCorrect - Checks whether the matrix A is suitable for iterative methods
%
% This function checks whether a given square matrix A is suitable for
% iterative methods such as the Jacobi or Gauss-Seidel method.
%
% The function first checks whether A is square. Then it computes the
% maximum of the normalized row sums and column sums (excluding the diagonal).
% If either maximum is less than 1, the matrix is considered suitable.
%
% If not, it checks the eigenvalues of A. If all eigenvalues have absolute
% value less than 1, the matrix is still considered suitable.
%
% Input:
%   A - a square matrix
%
% Output:
%   ret - logical value (true if A is acceptable, false otherwise)

Len = size(A,1);

if Len ~= size(A,2)
    disp("A is not square matrix!");
    ret = false;
    return;
end

MaxSum = 0;
MaxSumForColums = 0;

for i = 1:Len
    rowSum = 0;
    colSum = 0;
    for j = 1:Len
        if i == j
            continue
        end
        % Row-wise sum (excluding diagonal)
        rowSum = rowSum + abs(A(i,j)/A(i,i));
        % Column-wise sum (excluding diagonal)
        colSum = colSum + abs(A(j,i)/A(j,j));
    end
    % Update maximum values
    if MaxSum < rowSum
        MaxSum = rowSum;
    end
    if MaxSumForColums < colSum
        MaxSumForColums = colSum;
    end
end

% If either the maximum row or column sum is < 1, matrix is acceptable
if MaxSumForColums < 1 || MaxSum < 1
    ret = true;
    return;
else
    % Fallback: check eigenvalues
    eigenvalues = eig(A);
    disp("Checking the eigenvalues");
    for i = 1:length(eigenvalues)
        if abs(eigenvalues(i)) >= 1
            ret = false;
            return;
        end
    end
end

ret = true;

end
