function A = vandermonde(t, n)
t = t(:);
m = length(t);
A = zeros(m, n+1);
for j = 0:n
    A(:, j+1) = t.^j;
end
end
