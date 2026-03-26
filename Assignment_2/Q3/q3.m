clc; clear; close all;

%% (a) Gram-Schmidt orthonormality check
A_test = rand(6, 4);
[Q_test, ~] = myQR(A_test);
QtQ = Q_test' * Q_test;
fprintf('(a) Q^T Q (should be ~I):\n');
disp(round(QtQ, 10));

%% (b) Verify A = QR reconstruction
[Q_test, R_test] = myQR(A_test);
err = norm(A_test - Q_test * R_test, 'fro');
fprintf('(b) Reconstruction error ||A - QR||_F = %.2e\n\n', err);

%% (c) Timing comparison on 50x10 matrix
A = rand(50, 10);
t_mine    = timeit(@() myQR(A));
t_builtin = timeit(@() qr(A));
fprintf('(c) myQR runtime:      %.4f ms\n', t_mine * 1000);
fprintf('    Built-in qr runtime: %.4f ms\n\n', t_builtin * 1000);

figure;
bar([t_mine, t_builtin] * 1000);
set(gca, 'XTickLabel', {'myQR', 'Built-in qr'});
ylabel('Time (ms)'); title('QR Decomposition Runtime Comparison');
grid on;
saveas(gcf, 'timing.png');
close all;

%% (d) Least squares Ax = b
rng(42);
b = rand(50, 1);
[Q50, R50] = myQR(A);

x_mine   = R50 \ (Q50' * b);
x_matlab = A \ b;

res_mine   = norm(b - A * x_mine);
res_matlab = norm(b - A * x_matlab);
fprintf('(d) Residual norm (myQR):      %.6f\n', res_mine);
fprintf('    Residual norm (MATLAB \\):  %.6f\n\n', res_matlab);

figure;
idx = 1:50;
plot(idx, b, 'bo', 'MarkerSize', 5); hold on;
plot(idx, A * x_mine, 'r-', 'LineWidth', 1.5);
xlabel('Index'); ylabel('Value');
title('Least Squares Fit: b vs Ax_{hat}');
legend('b (data)', 'Ax_{hat} (fit)');
grid on;
saveas(gcf, 'lsq_fit.png');
close all;

fprintf('Q3 done\n');
