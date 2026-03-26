function [Q, R] = myQR(A)
[m, n] = size(A);
Q = zeros(m, n);

for k = 1:n
    v = A(:, k);
    for j = 1:k-1
        v = v - dot(Q(:,j), A(:,k)) * Q(:,j);
    end
    Q(:, k) = v / norm(v);
end

R = Q' * A;
end
