clc; clear; close all;

t = [0 1 2 3 4 5]';
y = [2.1 3.9 6.2 7.8 10.1 12.3]';

t_fine = linspace(0, 5, 200)';
norms = zeros(1, 2);
fits  = zeros(length(t_fine), 2);

for deg = [1 2]
    A = vandermonde(t, deg);

    c = inv(A' * A) * (A' * y);

    b_hat = A * c;
    res   = y - b_hat;

    fprintf('Degree %d fit coefficients: ', deg);
    fprintf('%.4f ', c); fprintf('\n');
    fprintf('  Residual norm: %.4f\n', norm(res));

    AtR = A' * res;
    fprintf('  A^T r (should be ~0): max abs = %.2e\n\n', max(abs(AtR)));

    norms(deg) = norm(res);
    A_fine = vandermonde(t_fine, deg);
    fits(:, deg) = A_fine * c;

    figure;
    plot(t, y, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); hold on;
    plot(t_fine, fits(:,deg), 'r-', 'LineWidth', 2);
    xlabel('t'); ylabel('y');
    title(sprintf('Polynomial Fit - Degree %d', deg));
    legend('Data', sprintf('Degree %d fit', deg));
    grid on;
    saveas(gcf, sprintf('fit_deg%d.png', deg));
    close;

    figure;
    stem(1:length(t), res, 'filled', 'LineWidth', 1.5);
    xlabel('Index'); ylabel('Residual');
    title(sprintf('Residual Vector - Degree %d', deg));
    grid on;
    saveas(gcf, sprintf('residual_deg%d.png', deg));
    close;
end

fprintf('Residual norms -> Degree 1: %.4f | Degree 2: %.4f\n', norms(1), norms(2));
if norms(2) < norms(1)
    fprintf('Degree 2 fits better (smaller residual norm).\n');
else
    fprintf('Degree 1 fits better (smaller residual norm).\n');
end

figure;
plot(t, y, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); hold on;
plot(t_fine, fits(:,1), 'r-', 'LineWidth', 2);
plot(t_fine, fits(:,2), 'g--', 'LineWidth', 2);
xlabel('t'); ylabel('y');
title('Degree 1 vs Degree 2 Polynomial Fit');
legend('Data', 'Degree 1', 'Degree 2');
grid on;
saveas(gcf, 'fit_comparison.png');
close;

fprintf('Q4 done\n');
